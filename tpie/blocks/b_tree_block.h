// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; c-file-style: "stroustrup"; -*-
// vi:set ts=4 sts=4 sw=4 noet :
// Copyright 2013, 2014, The TPIE development team
//
// This file is part of TPIE.
//
// TPIE is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// TPIE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TPIE.  If not, see <http://www.gnu.org/licenses/>

///////////////////////////////////////////////////////////////////////////////
/// \file blocks/b_tree_block.h  B+ tree internal node buffer view
///////////////////////////////////////////////////////////////////////////////

#ifndef TPIE_BLOCKS_B_TREE_BLOCK_H
#define TPIE_BLOCKS_B_TREE_BLOCK_H

#include <tpie/blocks/b_tree_bits.h>
#include <tpie/blocks/b_tree_leaf.h>

namespace tpie {

namespace blocks {

template <typename Key, typename Value, typename Compare,typename KeyExtract, typename Augment, typename Augmentor>
class b_tree_block {
public:
	static memory_size_type calculate_fanout(memory_size_type blockSize) {
		memory_size_type perChild = sizeof(block_handle) + sizeof(Augment);
		memory_size_type perKey =  sizeof(Key) + perChild;
		blockSize -= sizeof(b_tree_header);
		blockSize -= perChild; // one more child pointer than keys
		return blockSize / perKey; // floored division
	}

	b_tree_block(block_buffer & buffer, const b_tree_parameters & params)
		: m_params(params)
		, m_keyExtract()
		, m_augmentor()
	{
		char * children = buffer.get() + sizeof(b_tree_header);
		char * augments = children + params.nodeMax * sizeof(block_handle);
		char * keys = augments + params.nodeMax * sizeof(Augment);

		m_header = reinterpret_cast<b_tree_header *>(buffer.get());
		m_children = reinterpret_cast<block_handle *>(children);
		m_augments = reinterpret_cast<Augment *>(augments);
		m_keys = reinterpret_cast<Key *>(keys);
	}

	// Called by b_tree::insert after splitting the root into `left` and `right`.
	void new_root(const Key & k, block_handle left, const Augment & leftAugment, block_handle right, const Augment & rightAugment) {
		m_header->degree = 2;
		m_keys[0] = k;

		m_children[0] = left;
		m_augments[0] = leftAugment;

		m_children[1] = right;
		m_augments[1] = rightAugment;
	}

	void clear() {
		m_header->degree = 0;
	}

	// Internal helper used by b_tree_builder.
	void push_first_child(block_handle block, const Augment & augment) {
		if (!empty())
			throw exception("push_first_child: !empty");

		m_children[0] = block;
		m_augments[0] = augment;
		m_header->degree = 1;
	}

	// Internal helper used by b_tree_builder.
	void push_child(Key k, block_handle block, const Augment & augment) {
		if (full())
			throw exception("push_child: full");

		++m_header->degree;
		m_keys[m_header->degree - 2] = k;
		m_children[m_header->degree - 1] = block;
		m_augments[m_header->degree - 1] = augment;
	}

	memory_size_type degree() const {
		return static_cast<memory_size_type>(m_header->degree);
	}

	memory_size_type keys() const {
		return degree() - 1;
	}

	// calculates the augment for the given block
	Augment augment() const {
		return m_augmentor(m_augments, m_augments + degree());
	}

	// Definition 1, second bullet:
	// Except for the root, all nodes have degree
	// between nodeMin and nodeMax
	// (contain between nodeMin - 1 and nodeMax - 1 elements)
	bool full() const {
		return degree() == m_params.nodeMax;
	}

	bool underfull() const {
		return degree() < m_params.nodeMin;
	}

	bool empty() const {
		return degree() == 0;
	}

	Key key(memory_size_type idx) const {
		if (idx > keys()) throw exception("Block key: Index out of bounds");
		return m_keys[idx];
	}

	block_handle child(memory_size_type idx) const {
		if (idx > degree()) throw exception("Block child: Index out of bounds");
		return m_children[idx];
	}

	const Augment & augment(memory_size_type idx) const  {
		if (idx > degree()) throw exception("Block augment: Index out of bounds");
		return m_augments[idx];
	}

	// set the augment for the child with index idx
	void set_augment(memory_size_type idx, const Augment & augment) {
		if(idx > degree()) throw exception("Block augment: Index out of bounds");
		m_augments[idx] = augment;
	}

	// Called by b_tree::insert
	// Pre-condition: !full()
	void insert(memory_size_type i, Key k, block_handle leftChild, const Augment & leftAugment, block_handle rightChild, const Augment & rightAugment) {
		if (full()) throw exception("Insert on full block");

		m_children[i] = leftChild;
		m_augments[i] = leftAugment;

		block_handle handle = rightChild;
		Augment augment = rightAugment;
		while (i < keys()) {
			std::swap(m_children[i+1], handle);
			std::swap(m_augments[i+1], augment);
			std::swap(m_keys[i], k);
			++i;
		}
		m_children[i+1] = handle;
		m_augments[i+1] = augment;
		m_keys[i] = k;
		++m_header->degree;
	}

	// Called by b_tree::insert
	// Pre-condition: full()
	Key split_insert(memory_size_type insertIndex,
					 Key insertKey,
					 block_handle leftChild,
					 const Augment & leftAugment,
					 block_handle rightChild,
					 const Augment & rightAugment,
					 block_buffer & leftBuf,
					 block_buffer & rightBuf)
	{
		if (!full()) throw exception("split_insert on non-full block");

		typedef key_ops<Key> O;
		typedef typename O::ptr_type KeyPtr;

		std::vector<block_handle> children(degree()+1);
		std::vector<Augment> augments(degree()+1);
		std::vector<KeyPtr> keys(this->keys() + 1);

		{
			for (memory_size_type i = 0; i < this->keys(); ++i) {
				memory_size_type dest = i;
				if (insertIndex <= i) ++dest;
				children[dest] = m_children[i];
				augments[dest] = m_augments[i];
				keys[dest] = O::get_ptr(m_keys[i]);
			}
			children[degree()] = m_children[degree()-1];
			augments[degree()] = m_augments[degree()-1];

			keys[insertIndex] = O::get_ptr(insertKey);
			children[insertIndex] = leftChild;
			augments[insertIndex] = leftAugment;
			children[insertIndex+1] = rightChild;
			augments[insertIndex+1] = rightAugment;
		}

		Key midKey;

		{
			b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> left(leftBuf, m_params);
			b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> right(rightBuf, m_params);

			memory_size_type in = 0;
			memory_size_type out;
			for (out = 0; in*2 < keys.size(); ++out) {
				left.m_children[out] = children[in];
				left.m_augments[out] = augments[in];
				left.m_keys[out] = O::get_val(keys[in]);
				++in;
			}
			left.m_children[out] = children[in];
			left.m_augments[out] = augments[in];
			left.m_header->degree = static_cast<uint64_t>(out + 1);

			midKey = O::get_val(keys[in]);
			++in;

			for (out = 0; in < keys.size(); ++out) {
				right.m_children[out] = children[in];
				right.m_augments[out] = augments[in];
				right.m_keys[out] = O::get_val(keys[in]);
				++in;
			}
			right.m_children[out] = children[in];
			right.m_augments[out] = augments[in];
			right.m_header->degree = static_cast<uint64_t>(out + 1);
		}

		m_header->degree = 0;
		return midKey;
	}

	// Called by b_tree::erase
	// Returns fuse_merge or fuse_share.
	fuse_result fuse_leaves(memory_size_type rightIndex,
							block_buffer & leftBuf,
							block_buffer & rightBuf,
							const Compare & comp)
	{
		b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> left(leftBuf, m_params);
		b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> right(rightBuf, m_params);
		Key k;
		switch (left.fuse_with(right, k, comp)) {
			case fuse_merge:
				std::copy(m_keys + rightIndex,
						  m_keys + keys(),
						  m_keys + (rightIndex - 1));
				std::copy(m_children + (rightIndex + 1),
						  m_children + degree(),
						  m_children + rightIndex);
				std::copy(m_augments + (rightIndex + 1),
						  m_augments + degree(),
						  m_augments + rightIndex);
				--m_header->degree;
				return fuse_merge;
			case fuse_share:
				m_keys[rightIndex-1] = k;
				return fuse_share;
			default:
				throw exception("Unreachable statement");
		}
	}

	// Called by b_tree::erase
	// Returns fuse_merge or fuse_share.
	fuse_result fuse(memory_size_type rightIndex,
					 block_buffer & leftBuf,
					 block_buffer & rightBuf)
	{
		b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> left(leftBuf, m_params);
		b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> right(rightBuf, m_params);

		std::vector<Key> keys(left.keys() + 1 + right.keys());
		std::vector<block_handle> children(left.degree() + right.degree());
		std::vector<Augment> augments(left.degree() + right.degree());

		{
			memory_size_type output = 0;
			for (memory_size_type i = 0; i < left.keys(); ++i) {
				keys[output] = left.key(i);
				children[output] = left.child(i);
				augments[output] = left.augment(i);
				++output;
			}
			keys[output] = key(rightIndex-1);
			children[output] = left.child(left.keys());
			augments[output] = left.augment(left.keys());
			++output;
			for (memory_size_type i = 0; i < right.keys(); ++i) {
				keys[output] = right.key(i);
				children[output] = right.child(i);
				augments[output] = right.augment(i);
				++output;
			}
			children[output] = right.child(right.keys());
			augments[output] = right.augment(right.keys());
			++output;
		}

		if (children.size() <= m_params.nodeMax) {
			std::copy(keys.begin(),
					  keys.end(),
					  &left.m_keys[0]);
			std::copy(children.begin(),
					  children.end(),
					  &left.m_children[0]);
			std::copy(augments.begin(),
					  augments.end(),
					  &left.m_augments[0]);
			left.m_header->degree = static_cast<uint64_t>(children.size());

			std::copy(&m_keys[rightIndex],
					  &m_keys[this->keys()],
					  &m_keys[rightIndex-1]);
			std::copy(&m_children[rightIndex+1],
					  &m_children[degree()],
					  &m_children[rightIndex]);
			std::copy(&m_augments[rightIndex+1],
					  &m_augments[degree()],
					  &m_augments[rightIndex]);
			--m_header->degree;

			return fuse_merge;

		} else {

			memory_size_type half = children.size()/2;
			std::copy(keys.begin(),
					  keys.begin() + (half - 1),
					  &left.m_keys[0]);
			std::copy(children.begin(),
					  children.begin() + half,
					  &left.m_children[0]);
			std::copy(augments.begin(),
					  augments.begin() + half,
					  &left.m_augments[0]);
			left.m_header->degree =
				static_cast<uint64_t>(half);

			m_keys[rightIndex-1] = keys[half-1];

			std::copy(keys.begin() + half,
					  keys.end(),
					  &right.m_keys[0]);
			std::copy(children.begin() + half,
					  children.end(),
					  &right.m_children[0]);
			std::copy(augments.begin() + half,
					  augments.end(),
					  &right.m_augments[0]);
			right.m_header->degree =
				static_cast<uint64_t>(children.size() - half);

			return fuse_share;
		}
	}

private:
	b_tree_header * m_header;
	block_handle * m_children;
	Augment * m_augments;
	Key * m_keys;
	b_tree_parameters m_params;
	KeyExtract m_keyExtract;
	Augmentor m_augmentor;
};

} // namespace blocks

} // namespace tpie

#endif // TPIE_BLOCKS_B_TREE_BLOCK_H
