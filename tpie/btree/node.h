// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :
// Copyright 2014 The TPIE development team
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

#ifndef _TPIE_BTREE_NODE_H_
#define _TPIE_BTREE_NODE_H_

#include <tpie/portability.h>
#include <tpie/tpie_assert.h>
#include <tpie/btree/base.h>
#include <boost/iterator/iterator_facade.hpp>
#include <vector>

namespace tpie {

/**
 * Type that is usefull for navigating a btree
 *
 * S is the type of the store used
 */
template <typename S>
class btree_node {
public:
	/**
	 * Type of the key of a value
	 */
	typedef typename S::key_type key_type;

	/**
	 * Type of the augment of a set of notes/values
	 */
	typedef typename S::augment_type augment_type;

	/**
	 * Type of values
	 */
	typedef typename S::value_type value_type;

	/**
	 * \brief Check if this note has a parent
	 *
	 * True iff this is not the root
	 */
	bool has_parent() const {
		if (m_is_leaf)
			return !m_path.empty();
		else
			return m_path.size() > 1;
	}

	/**
	 * \brief Move to the parent node
	 *
	 * Requires hasParent()
	 */
	void parent() {
		if (m_is_leaf)
			m_is_leaf = false;
		else
			m_path.pop_back();
	}
		
	/**
	 * \brief Move to the ith child
	 *
	 * Requires !leaf() and i < count()
	 */
	void child(size_t i) {
		tp_assert(!m_is_leaf, "Is leaf");
		tp_assert(i < count(), "Invalid i");
		if (m_path.size() + 1 == m_store->height()) {
			m_is_leaf = true;
			m_leaf = m_store->get_child_leaf(m_path.back(), i);
		} else
			m_path.push_back(m_store->get_child_internal(m_path.back(), i));
	}
	
	/**
	 * \brief Return the parent node
	 *
	 * Requires hasParent()
	 */
	btree_node get_parent() const {
		btree_node n=*this;
		n.parent();
		return n;
	}
	
	/**
	 * \brief Return the ith child note
	 *
	 * Requires !leaf() and i < count()
	 */
	btree_node get_child(size_t i) const {
		btree_node n=*this;
		n.child(i);
		return n;
	}
	
	/**
	 * \brief Return true if this is a leaf note
	 */
	bool leaf() const {
		return m_is_leaf;
	}
	
	/**
	 * \brief Return the augment of the ith child
	 *
	 * Requires !leaf()
	 */
	const augment_type & augment(size_t i) const {
		return m_store->augment(m_path.back(), i);
	}
	
	/**
	 * \brief Return the minimal key of the i'th child
	 */
	key_type min_key(size_t i) const {
		if (m_is_leaf)
			return m_store->min_key(m_leaf, i);
		else
			return m_store->min_key(m_path.back(), i);
	}

	/**
	 * \brief Return the i'th value
	 *
	 * Requires leaf()
	 */
	const value_type & value(size_t i) const {
		tp_assert(m_is_leaf, "Not leaf");
		return m_store->get(m_leaf, i);
	}
	
	/**
	 * \brief Return the i'th value
	 *
	 * Requires leaf()
	 */
	value_type & value(size_t i) {
		tp_assert(m_is_leaf, "Not leaf");
		return m_store->get(m_leaf, i);
	}
	
	/**
	 * \brief Return the number of children or values
	 */
	size_t count() const {
		if (m_is_leaf)
			return m_store->count(m_leaf);
		else
			return m_store->count(m_path.back());
	}

	/**
	 * \brief Return the index of this node in its parent
	 *
	 * Requires has_parent()
	 */
	size_t index() const {
		if (m_is_leaf)
			return m_store->index(m_leaf, m_path.back());
		return m_store->index(m_path.back(), m_path[m_path.size()-2]);
	}
	
	btree_node(): m_store(NULL) {}
private:
	typedef typename S::leaf_type leaf_type;
	typedef typename S::internal_type internal_type;

	btree_node(S * store, leaf_type root)
		: m_store(store), m_leaf(root), m_is_leaf(true) {
	}

	btree_node(S * store, internal_type root)
		: m_store(store), m_is_leaf(false) {
		m_path.push_back(root);
	}

	btree_node(S * store, std::vector<internal_type> path, leaf_type leaf)
		: m_store(store), m_path(path), m_leaf(leaf), m_is_leaf(true) {
	}

	S * m_store;
	std::vector<internal_type> m_path;
	leaf_type m_leaf;
	bool m_is_leaf;

	template <typename, typename, typename>
	friend class btree;

	template <typename>
	friend class btree_iterator;
};

template <typename S>
class btree_iterator: public boost::iterator_facade<
	btree_iterator<S>,
	typename S::value_type,
	boost::bidirectional_traversal_tag> {
private:
	typedef typename S::internal_type internal_type;
	typedef typename S::leaf_type leaf_type;
	typedef typename S::value_type value_type;
	typedef typename S::key_type key_type;

	S * m_store;
	std::vector<internal_type> m_path;
	size_t m_index;
	leaf_type m_leaf;

	template <typename, typename, typename>
	friend class btree;

	btree_iterator(S * store): m_store(store) {}

	void goto_item(const std::vector<internal_type> & p, leaf_type l, size_t i) {
		m_path = p;
		m_leaf = l;
		m_index = i;
	}


	void goto_begin() {
		m_path.clear();
		if (m_store->height() == 1) {
			m_leaf = m_store->get_root_leaf();
			m_index = 0;
			return;
		}
		internal_type n = m_store->get_root_internal();
		for (size_t i=2;; ++i) {
			m_path.push_back(n);
			if (i == m_store->height()) {
				m_leaf = m_store->get_child_leaf(n, 0);
				m_index = 0;
				return;
			}
			n = m_store->get_child_internal(n, 0);
		}		
	}

	void goto_end() {
		m_path.clear();
		if (m_store->height() == 1) {
			m_leaf = m_store->get_root_leaf();
			m_index = m_store->count(m_leaf);
			return;
		}
		internal_type n = m_store->get_root_internal();
		for (size_t i=2;; ++i) {
			m_path.push_back(n);
			if (i == m_store->height()) {
				m_leaf = m_store->get_child_leaf(n, m_store->count(n)-1);
				m_index = m_store->count(m_leaf);
				return;
			}
			n = m_store->get_child_internal(n, m_store->count(n)-1);
		}		
	}


public:
	btree_iterator(): m_index(0), m_leaf() {}

	value_type & dereference() const { 
		return m_store->get(m_leaf, m_index);
	}

	bool equal(const btree_iterator & o) const {
		return m_index == o.m_index && m_leaf == o.m_leaf;
	}
	
	size_t index() const {return m_index;}

	btree_node<S> leaf() const {
		return btree_node<S>(m_store, m_path, m_leaf);
	}
	
	void decrement() {
		if (m_index > 0) {
			--m_index;
			return;
		}
		
		size_t i=m_store->index(m_leaf, m_path.back());
		size_t x=0;
		while (i == 0) {
			internal_type n = m_path.back();
			m_path.pop_back();
			i = m_store->index(n, m_path.back());
			++x;
		}
		--i;

		while (x != 0) {
			m_path.push_back(m_store->get_child_internal(m_path.back(), i));
			i = m_store->count(m_path.back())-1;							 
			--x;
		}
		
		m_leaf = m_store->get_child_leaf(m_path.back(), i);
		m_index = m_store->count(m_leaf)-1;
	}

	void increment() {
		m_index++;
		if (m_index < m_store->count(m_leaf)) return;
		if (m_path.empty()) return; //We are at the end

		size_t i=m_store->index(m_leaf, m_path.back());
		size_t x=0;
		while (i +1 == m_store->count(m_path.back())) {
			internal_type n = m_path.back();
			m_path.pop_back();
			if (m_path.empty()) {
				goto_end();
				return;
			}
			i = m_store->index(n, m_path.back());
			++x;
		}
		++i;
		while (x != 0) {
			m_path.push_back(m_store->get_child_internal(m_path.back(), i));
			i = 0;
			--x;
		}
		m_leaf = m_store->get_child_leaf(m_path.back(), i);
		m_index = 0;
	}

};


} //namespace tpie
#endif //_TPIE_BTREE_NODE_H_
