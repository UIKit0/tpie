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
/// \file blocks/b_tree_bits.h  B+ tree helper classes.
///////////////////////////////////////////////////////////////////////////////

#ifndef TPIE_BLOCKS_B_TREE_BITS_H
#define TPIE_BLOCKS_B_TREE_BITS_H

#include <tpie/blocks/block_collection.h>

namespace tpie {

namespace blocks {

template<typename T>
struct identity_key_extract {
	T operator()(const T & t) const {
		return t;
	}
};

struct empty_augment {

};

template <typename T>
struct empty_augmentor {
	empty_augment operator()(empty_augment *, empty_augment *) { // calculate the augment for an internal node
		return empty_augment();
	};

	empty_augment operator()(T *, T *) { // calculate the augment for a leaf node
		return empty_augment();
	}
};

struct b_tree_header {
	/** Number of child blocks or leaf values. */
	uint64_t degree;
};

// Contains the "branching parameter" nodeMax and the "leaf parameter" leafMax,
// which are called b and k respectively in the paper.
struct b_tree_parameters {
	uint64_t nodeMin;
	uint64_t nodeMax;
	uint64_t leafMin;
	uint64_t leafMax;
};

// key_ops is used in b_tree_block::split_insert to decide whether
// to take pointers into the original key arrays (expensive due to indirection)
// or to copy keys into the temporary storage (expensive in space).
// If the key is as big as its pointer, it makes no sense to introduce indirection.
// On the other hand, if the key is larger than its pointer, we don't want
// to waste memory on copying keys around.
template <typename Key,
		 bool PtrSized = sizeof(Key) <= sizeof(Key *)>
class key_ops {
public:
	typedef Key * ptr_type;
	static ptr_type get_ptr(Key & k)    { return &k; }
	static Key      get_val(ptr_type p) { return *p; }
};

template <typename Key>
class key_ops<Key, true> {
public:
	typedef Key ptr_type;
	static ptr_type get_ptr(Key & k)    { return k; }
	static Key      get_val(ptr_type p) { return p; }
};

// Result of a fuse operation.
enum fuse_result {
	/** `left` and `right` are still in use. */
	fuse_share,

	/** All of `right` was merged into `left`. */
	fuse_merge
};

// Functor for partitioning an array of Values according to a given Key.
template <typename Key, typename Value, typename Compare, typename KeyExtract, typename Augment, typename Augmentor>
class key_less_than {
	Compare m_comp;
	KeyExtract m_keyExtract;
	Key m_key;

public:
	key_less_than(const Compare & comp, const Value & v)
		: m_comp(comp)
		, m_keyExtract()
		, m_key(m_keyExtract(v))
	{
	}

	bool operator()(const Value & v) const {
		return m_comp(m_keyExtract(v), m_key);
	}
};

// Comparator for sorting an array of Values by their Keys.
template <typename Key, typename Value, typename Compare, typename KeyExtract, typename Augment, typename Augmentor>
class key_less {
	Compare m_comp;
	KeyExtract m_keyExtract;

public:
	key_less(const Compare & comp)
		: m_comp(comp)
		, m_keyExtract()
	{
	}

	bool operator()(const Value & v1, const Value & v2) const {
		return m_comp(m_keyExtract(v1), m_keyExtract(v2));
	}
};

// A path in the B tree of height h (having the root at level 0 and leaves at level h)
// is a sequence ((b_0, i_0), (b_1, i_1), ..., (b_(h-1), i_(h-1)))
// where b_0 is the block handle of the root block,
// and b_(j+1) is the i_j'th block referred to by block b_j.
// A b_tree_path supports the stack operation PUSH through the `follow` method,
// and the stack operation POP through the `parent` method.
class b_tree_path {
public:
	// Pre-condition:
	// If empty(), b refers to the root of the B tree.
	// If !empty(), b is the index'th child of block `current_block()`.
	void follow(block_handle b, memory_size_type index) {
		m_components.push_back(std::make_pair(b, index));
	}

	// Pre-condition: !empty()
	void parent() {
		m_components.pop_back();
	}

	block_handle current_block() const {
		return m_components[m_components.size()-1].first;
	}

	memory_size_type current_index() const {
		return m_components[m_components.size()-1].second;
	}

	bool empty() const {
		return m_components.empty();
	}

private:
	std::vector<std::pair<block_handle, memory_size_type> > m_components;
};

// Finite state machine embedded in the b_tree_builder
// to check that methods are called in the right order.
struct builder_state {
	enum type {
		/** No items have been pushed. */
		EMPTY,
		/** Items have been pushed, but end() has not been called. */
		BUILDING,
		/** end() has been called. */
		BUILT
	};
};

} // namespace blocks

} // namespace tpie

#endif // TPIE_BLOCKS_B_TREE_BITS_H
