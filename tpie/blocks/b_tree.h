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
/// \file blocks/b_tree.h  B+ tree in external memory
///
/// The implementation is based off the B tree described in
/// Lars Arge, External Memory Geometric Data Structures (2005).
///
/// The main class `b_tree` handles reading from and writing to the underlying
/// block collection.
///
/// The helper classes `b_tree_block` and `b_tree_leaf` operate on individual
/// blocks that have been read and support primitive operations such as
/// inserting, erasing, splitting and fusing.
///
/// The `b_tree` class uses these to insert and erase elements on a tree-wide
/// scale.
///
/// There is a clear division of responsibility: The block and leaf classes
/// have no access to reading from and writing to the block collection,
/// and the b_tree class cannot modify the key/child/value arrays directly,
/// but must do so through the helpers.
///////////////////////////////////////////////////////////////////////////////

#ifndef TPIE_BLOCKS_B_TREE_H
#define TPIE_BLOCKS_B_TREE_H

#include <boost/filesystem.hpp>

#include <tpie/tempname.h>
#include <tpie/blocks/block_collection.h>
#include <tpie/blocks/b_tree_bits.h>
#include <tpie/blocks/b_tree_leaf.h>
#include <tpie/blocks/b_tree_block.h>

namespace tpie {

namespace blocks {

template <typename Key, typename Value, typename Compare, typename KeyExtract, typename Augment, typename Augmentor>
class b_tree_builder;

template <typename Key, typename Value=Key, typename Compare=std::less<Key>, typename KeyExtract=identity_key_extract<Value>, typename Augment=empty_augment, typename Augmentor=empty_augmentor<Value> >
class b_tree {
public:
	b_tree()
		: m_root(0)
		, m_treeHeight(0)
		, m_keyExtract()
		, m_augmentor()
	{
		set_default_parameters();
	}

	~b_tree() {
		close();
	}

	bool is_open() const {
		return m_blocks.is_open();
	}

	void set_default_parameters() {
		if (is_open()) throw exception("set_default_parameters: block collection already open");
		m_params.nodeMax = b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor>::calculate_fanout(block_size());
		m_params.nodeMin = (m_params.nodeMax + 3)/4;
		m_params.leafMax = b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor>::calculate_fanout(block_size());
		m_params.leafMin = (m_params.leafMax + 3)/4;
		verify_parameters();
		log_parameters();
	}

	b_tree_parameters get_parameters() const {
		return m_params;
	}

	void set_parameters(const b_tree_parameters & params) {
		if (is_open()) throw exception("set_parameters: block collection already open");
		b_tree_parameters prev = m_params;
		m_params = params;
		try {
			verify_parameters();
			log_parameters();
		} catch (const exception &) {
			m_params = prev;
			throw;
		}
	}

private:
	void verify_parameters() {
		if (m_params.nodeMin < 2)
			throw exception("Block size too small; nodeMin >= 2 violated");
		if (m_params.nodeMax < m_params.nodeMin*2-1)
			throw exception("Block size too small; nodeMax >= 2a-1 violated");
		if (m_params.leafMin < 2)
			throw exception("Block size too small; leafMin >= 2 violated");
		if (m_params.leafMax < m_params.leafMin*2-1)
			throw exception("Block size too small; leafMax >= 2a-1 violated");
	}

	void log_parameters() {
		log_debug()
			<< "B tree parameters\n"
			<< "Node degree in [" << m_params.nodeMin << ", " << m_params.nodeMax << "]\n"
			<< "Leaf degree in [" << m_params.leafMin << ", " << m_params.leafMax << "]\n"
			<< std::flush;
	}

	void augment_path(b_tree_path p, Augment a) { // updates augments for all nodes in the given path
		block_buffer buf;
		while(!p.empty()) {
			m_blocks.read_block(p.current_block(), buf);
			b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> block(buf, m_params);
			block.set_augment(p.current_index(), a);
			m_blocks.write_block(buf);

			a = block.augment();
			p.parent();
		}
	}
public:
	///////////////////////////////////////////////////////////////////////////
	/// \brief  Insert value into the B tree.
	///////////////////////////////////////////////////////////////////////////
	void insert(Value v) {
		if (!is_open()) throw exception("insert: block collection not open");
		block_buffer buf;
		Key k = m_keyExtract(v);
		// Find the leaf in which the value should be inserted.
		b_tree_path p = key_path(buf, k);
		block_handle leftChild;
		block_handle rightChild;

		Augment leftAugment; // the augment for leftChild
		Augment rightAugment; // the augment for rightChild

		{
			b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> leaf(buf, m_params);

			// If the leaf is not full, we do a cheap insert, augment the path and return.
			if (!leaf.full()) {
				leaf.insert(v);
				m_blocks.write_block(buf);
				augment_path(p, leaf.augment());
				return;
			}

			// Split the leaf.
			block_buffer & left = buf;
			block_buffer right;

			m_blocks.get_free_block(right);

			k = leaf.split_insert(v, right, m_comp);
			m_blocks.write_block(left);
			m_blocks.write_block(right);

			// Proceed with recursive insertion below.
			b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> leftLeaf(left, m_params);
			leftChild = left.get_handle();
			leftAugment = leftLeaf.augment();

			b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> rightLeaf(right, m_params);
			rightChild = right.get_handle();
			rightAugment = rightLeaf.augment();
		}

		if (p.empty()) {
			// Special case: The root was previously a single leaf
			// which has now been split into two.
			// There is no need to augment the path since it is empty.
			m_blocks.get_free_block(buf);
			b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> block(buf, m_params);
			block.new_root(k, leftChild, leftAugment, rightChild, rightAugment);
			m_blocks.write_block(buf);
			++m_treeHeight;
			m_root = buf.get_handle();
			log_debug() << "Increase tree height to " << m_treeHeight << "; root is now " << m_root << std::endl;
			return;
		}

		b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> block(buf, m_params);

		// Repeatedly split blocks until we hit a non-full block
		// or we hit the root of the tree.
		while (!p.empty()) {
			m_blocks.read_block(p.current_block(), buf);
			if (!block.full()) break;
			// else, we split the block
			block_buffer left;
			block_buffer right;
			m_blocks.get_free_block(left);
			m_blocks.get_free_block(right);
			k = block.split_insert(p.current_index(), k, leftChild, leftAugment, rightChild, rightAugment, left, right);
			m_blocks.write_block(left);
			m_blocks.write_block(right);

			b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> leftLeaf(left, m_params);
			leftChild = left.get_handle();
			leftAugment = leftLeaf.augment();

			b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> rightLeaf(left, m_params);
			rightChild = right.get_handle();
			rightAugment = rightLeaf.augment();

			p.parent();
			m_blocks.free_block(buf);
		}

		if (p.empty()) {
			// We split the root and do not augment the path since it is empty.
			m_blocks.get_free_block(buf);
			block.new_root(k, leftChild, leftAugment, rightChild, rightAugment);
			m_blocks.write_block(buf);
			++m_treeHeight;
			m_root = buf.get_handle();
			log_debug() << "Increase tree height to " << m_treeHeight << "; root is now " << m_root << std::endl;
		} else {
			// We hit a non-full block.
			m_blocks.read_block(p.current_block(), buf);
			block.insert(p.current_index(), k, leftChild, leftAugment, rightChild, rightAugment);
			m_blocks.write_block(buf);
			p.parent();

			augment_path(p, block.augment());
		}
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Erase value from B tree given its key.
	///////////////////////////////////////////////////////////////////////////
	void erase(Key k) {
		if (!is_open()) throw exception("erase: block collection not open");
		block_buffer buf;
		// Find leaf from which to erase.
		b_tree_path p = key_path(buf, k);
		{
			b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> leaf(buf, m_params);

			leaf.erase(k, m_comp);

			m_blocks.write_block(buf);
			// If leaf is not underfull, or the leaf is the root of the tree,
			// we augment the path and quickly return.
			if (p.empty() || !leaf.underfull()) {
				augment_path(p, leaf.augment());
				return;
			}
		}

		// Recursively fuse blocks until we hit a block that is not underfull.

		memory_size_type rightIndex = p.current_index() == 0 ? 1 : p.current_index();
		block_buffer left;
		block_buffer right;

		m_blocks.read_block(p.current_block(), buf);
		b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> block(buf, m_params);

		m_blocks.read_block(block.child(rightIndex-1), left);
		m_blocks.read_block(block.child(rightIndex), right);
		switch (block.fuse_leaves(rightIndex, left, right, m_comp)) {
			case fuse_share:
				//log_debug() << "fuse_leaves() == fuse_share" << std::endl;
				m_blocks.write_block(buf);
				m_blocks.write_block(left);
				m_blocks.write_block(right);
				return;
			case fuse_merge:
				//log_debug() << "fuse_leaves() == fuse_merge" << std::endl;
				m_blocks.write_block(buf);
				m_blocks.write_block(left);
				m_blocks.free_block(right);
				break;
		}

		p.parent();
		Augment a =  block.augment();

		while (!p.empty() && block.underfull()) {
			memory_size_type i = p.current_index();
			m_blocks.read_block(p.current_block(), buf);
			block.set_augment(i, a);

			memory_size_type rightIndex = (i == 0) ? 1 : i;
			m_blocks.read_block(block.child(rightIndex-1), left);
			m_blocks.read_block(block.child(rightIndex), right);
			switch (block.fuse(rightIndex, left, right)) {
				case fuse_share:
					//log_debug() << "Erase fuse_share of " << left.get_handle()
						//<< " and " << right.get_handle() << std::endl;
					m_blocks.write_block(buf);
					m_blocks.write_block(left);
					m_blocks.write_block(right);
					a = block.augment();
					p.parent();
					augment_path(p, a);
					return;
				case fuse_merge:
					//log_debug() << "Erase fuse_merge of " << left.get_handle()
						//<< " and " << right.get_handle() << std::endl;
					m_blocks.write_block(buf);
					m_blocks.write_block(left);
					m_blocks.free_block(right);
					a = block.augment();
					break;
			}
			p.parent();
		}
		if (p.empty() && block.degree() == 1) {
			m_root = block.child(0);
			m_blocks.free_block(buf);
			--m_treeHeight;
			log_debug() << "Decrease tree height to " << m_treeHeight << "; root is now " << m_root << std::endl;
			return;
		}

		augment_path(p, a);
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Search B tree for item with given Key.
	///
	/// Returns 1 if found, and 0 if not found.
	///////////////////////////////////////////////////////////////////////////
	memory_size_type count(Key k) {
		if (!is_open()) throw exception("count: block collection not open");
		block_buffer buf;
		b_tree_path p = key_path(buf, k);
		b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> leaf(buf, m_params);
		return leaf.count(k, m_comp);
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Search B tree for item with given Key.
	///
	/// If found, assigns the result to *out. `out` must be non-null.
	/// If you do not need the value associated to the key,
	/// use count() instead.
	///
	/// Returns true if found; false otherwise.
	///////////////////////////////////////////////////////////////////////////
	bool try_find(Key k, Value * out) {
		if (!is_open()) throw exception("try_find: block collection not open");
		block_buffer buf;
		b_tree_path p = key_path(buf, k);
		b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> leaf(buf, m_params);
		memory_size_type i = leaf.index_of(k, m_comp);
		if (i == leaf.degree()) {
			return false;
		} else {
			*out = leaf[i];
			return true;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Find value associated to given key.
	///
	/// Precondition: count(k) == 1.
	///
	/// If you are not sure if the tree contains a value for the given key,
	/// use `try_find` instead.
	///////////////////////////////////////////////////////////////////////////
	Value find(Key k) {
		if (!is_open()) throw exception("find: block collection not open");
		Value out;
		if (try_find(k, &out))
			return out;
		else
			throw exception("Value with given key not found");
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Open an anonymous B tree.
	///////////////////////////////////////////////////////////////////////////
	void open() {
		open_inner(m_tempFile.path());
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Open the B tree to the given block collection.
	///////////////////////////////////////////////////////////////////////////
	void open(tpie::temp_file & tempFile) {
		open_inner(tempFile.path());
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Open the B tree to the given block collection.
	///////////////////////////////////////////////////////////////////////////
	void open(std::string path) {
		open_inner(path);
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Close B tree.
	///////////////////////////////////////////////////////////////////////////
	void close() {
		if (is_open()) {
			m_blocks.close();
			m_tempFile.free();
		}
	}

private:
	void open_inner(const std::string & path) {
		if (is_open()) throw exception("open: block collection already open");
		m_blocks.open(path, true);
	}

	memory_size_type block_size() {
		return block_collection::default_block_size();
	}

	void read_root(block_buffer & b) {
		if (m_root == block_handle(0)) {
			m_root = m_blocks.get_free_block();
			m_treeHeight = 0;
		}
		m_blocks.read_block(m_root, b);
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Search the tree for the insertion point of a given key.
	///
	/// The leaf in which to insert is read into `buf`, and the path to the
	/// buffer is returned.
	///////////////////////////////////////////////////////////////////////////
	b_tree_path key_path(block_buffer & buf, Key k) {
		b_tree_path res;

		read_root(buf);

		for (memory_size_type i = 0; i < m_treeHeight; ++i) {
			b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> b(buf, m_params);

			memory_size_type j;
			for (j = 0; j != b.keys(); ++j)
				if (m_comp(k, b.key(j))) break;

			res.follow(buf.get_handle(), j);

			if (b.child(j) == block_handle(0))
				throw exception("Child pointer is 0 in non-leaf");
			else
				m_blocks.read_block(b.child(j), buf);
		}

		return res;
	}

public:
	struct leaf_range {
		Value const * begin;
		Value const * end;
	};

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Iterate through values of the tree in the range [a;b].
	///
	/// \tparam It  A functor that accepts a leaf_range as a param
	///////////////////////////////////////////////////////////////////////////
	template <typename Func>
	void leaf_range_report(const Key & a, const Key & b, Func func) {
		if(!is_open()) throw exception("in_order_dump: block collection not open");
		if(m_root == block_handle(0)) {
			log_debug() << "in_order_dump: Empty tree" << std::endl;
			return;
		}
		leaf_range_report_visit(a, b, func, m_root, m_treeHeight);
	}

private:
	template <typename Func>
	void leaf_range_report_visit(const Key & a, const Key & b, Func func, block_handle id, memory_size_type leafDistance) {
		if(id == block_handle(0)) return;
		block_buffer buf;
		m_blocks.read_block(id, buf);
		if(leafDistance == 0) { // the node is a leaf
			b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> leaf(buf, m_params);

			leaf_range args;
			args.begin = leaf.begin();
			args.end  = leaf.end();
			func(args);
			return;
		}

		// the node is an internal node
		b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> block(buf, m_params);
		if(block.underfull() && id != m_root)
			log_error() << "in_order_dump: Underfull non-root block " << id << std::endl;

		memory_size_type i, j;
		// find first key k where a <= k
		for(i = 0; i != block.keys(); ++i)
			if(!m_comp(a, block.key(i))) break;

		// find the first key k where !(k <= b)
		for(j = i; j != block.keys(); ++j)
			if(m_comp(b, block.key(j))) break;


		for(memory_size_type l = i; l <= j; ++l) {
			leaf_range_report_visit(a, b, func, block.child(l), leafDistance-1);
		}
	}

	template <typename It>
	class it_functor {
	public:
		it_functor(It it, const Key & a, const Key & b)
		: it(it)
		, m_comp()
		, a(a)
		, b(b)
		{}

		void operator()(leaf_range range) {
			while(range.begin != range.end) {
				if(!m_comp(*range.begin, a) && !m_comp(b, *range.begin))
					*it = *(range.begin);
				++it;
				++(range.begin);
			}
		}
	private:
		It it;
		Compare m_comp;
		const Key & a;
		const Key & b;
	};

public:
	///////////////////////////////////////////////////////////////////////////
	/// \brief  Iterate through values of the tree in the range [a;b].
	///
	/// \tparam it  An output iterator type that has Value among its
	///             value types.
	///////////////////////////////////////////////////////////////////////////
	template <typename It>
	void range_report(const Key & a, const Key & b, It it) {
		it_functor<It> func(it, a, b);
		leaf_range_report(a, b, func);
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Returns the augment for the values in the interval [a;b]
	///////////////////////////////////////////////////////////////////////////////
	Augment augment(const Key & a, const Key & b) {
		if(!is_open()) throw exception("in_order_dump: block collection not open");
		if(m_root == block_handle(0)) {
			log_debug() << "in_order_dump: Empty tree" << std::endl;
			return Augment();
		}
		return augment_visit(a, b, false, false, m_root, m_treeHeight);
	}
private:
	Augment augment_visit(const Key & a, const Key & b, bool leftInside, bool rightInside, block_handle id, memory_size_type leafDistance) {
		// leftInside is true if a <= e for all keys in the tree
		// rightInside is true if e <= b for all keys in the tree

		if(id == block_handle(0)) {
			log_debug() << "augment_visit: discovered block with block_handle 0" << std::endl;
			return Augment();
		}
		block_buffer buf;
		m_blocks.read_block(id, buf);
		if(leafDistance == 0) { // the node is a leaf
			b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> leaf(buf, m_params);

			std::vector<Value> values;
			for(Value const * i = leaf.begin(); i != leaf.end(); ++i) {
				if(!m_comp(*i, a) && !m_comp(b, *i))
					values.push_back(*i);
			}

			Value * begin = &values[0];
			return m_augmentor(begin, begin + values.size());
		}

		// the node is an internal node
		b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> block(buf, m_params);
		if(block.underfull() && id != m_root)
			log_error() << "in_order_dump: Underfull non-root block " << id << std::endl;

		memory_size_type i, j;
		// find first key k where a <= k
		for(i = 0; i != block.keys(); ++i)
			if(!m_comp(block.key(i), a)) break;

		// find the first key k where b < k
		for(j = i; j != block.keys(); ++j)
			if(m_comp(b, block.key(j))) break;


		Augment augments[j-i+1];
		for(memory_size_type l = i; l <= j; ++l) {
			bool newLeftInside = leftInside;
			newLeftInside |= (l > 0 && !m_comp(block.key(l-1), a)); // inside if the previous key k satisfy a <= k
			bool newRightInside = rightInside;
			newRightInside |= (l < (block.degree() - 1) && !m_comp(b, block.key(l))); // inside if the next key k satisfy k <= b
			if(newLeftInside && newRightInside)
				augments[l-i] = block.augment(l);
			else
				augments[l-i] = augment_visit(a, b, newLeftInside, newRightInside, block.child(l), leafDistance-1);
		}

		return m_augmentor(augments, augments+(j-i+1));
	}
public:

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Iterate through values of the tree in-order.
	///
	/// \tparam It  An output iterator type that has Value among its
	///             value types.
	///////////////////////////////////////////////////////////////////////////
	template <typename It>
	void in_order_dump(It it) {
		if (!is_open()) throw exception("in_order_dump: block collection not open");
		if (m_root == block_handle(0)) {
			log_debug() << "in_order_dump: Empty tree" << std::endl;
			return;
		}
		in_order_dump_visit(it, m_root, m_treeHeight);
	}

private:
	template <typename It>
	void in_order_dump_visit(It it, block_handle id, memory_size_type leafDistance) {
		if (id == block_handle(0)) return;
		block_buffer buf;
		m_blocks.read_block(id, buf);
		if (leafDistance == 0) {
			b_tree_leaf<Key, Value, Compare, KeyExtract, Augment, Augmentor> leaf(buf, m_params);
			std::vector<Value> vals;
			vals.reserve(leaf.degree());
			for (memory_size_type i = 0; i < leaf.degree(); ++i) {
				vals.push_back(leaf[i]);
			}
			std::sort(vals.begin(), vals.end(), key_less<Key, Value, Compare, KeyExtract, Augment, Augmentor>(m_comp));
			for (memory_size_type i = 0; i < vals.size(); ++i) {
				*it = vals[i];
				++it;
			}
		} else {
			b_tree_block<Key, Value, Compare, KeyExtract, Augment, Augmentor> block(buf, m_params);
			if (block.underfull() && id != m_root) {
				log_error() << "in_order_dump: Underfull non-root block " << id << std::endl;
			}
			for (memory_size_type i = 0; i < block.degree(); ++i) {
				in_order_dump_visit(it, block.child(i), leafDistance-1);
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Internal b_tree_builder helper.
	///////////////////////////////////////////////////////////////////////////
	void set_root(block_handle root, memory_size_type treeHeight) {
		m_root = root;
		m_treeHeight = treeHeight;
	}

	tpie::temp_file m_tempFile;
	block_collection m_blocks;
	block_handle m_root;
	Compare m_comp;
	memory_size_type m_treeHeight;
	b_tree_parameters m_params;
	KeyExtract m_keyExtract;
	Augmentor m_augmentor;

	friend class b_tree_builder<Key, Value, Compare, KeyExtract, Augment, Augmentor>;
};

} // namespace blocks

} // namespace tpie

#include <tpie/blocks/b_tree_builder.h>

#endif // TPIE_BLOCKS_B_TREE_H
