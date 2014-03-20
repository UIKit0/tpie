// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; c-file-style: "stroustrup"; -*-
// vi:set ts=4 sts=4 sw=4 noet cino+=(0 :
// Copyright 2014, The TPIE development team
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
#ifndef __TPIE_INTERNAL_DEQUE_H__
#define __TPIE_INTERNAL_DEQUE_H__

///////////////////////////////////////////////////////////////////////////////
/// \file internal_deque.h
/// Generic internal deque with known memory requirements.
///////////////////////////////////////////////////////////////////////////////
#include <tpie/array.h>
#include <tpie/util.h>
#include <tpie/tpie_assert.h>

namespace tpie {

///////////////////////////////////////////////////////////////////////////////
/// \brief A generic internal circular deque
///
/// The deque supports a fixed number of unpopped elements between
/// calls to clear. The number of elements is given as an argument
/// to the constructor or to resize.
///
/// \tparam T The type of items stored in the deque
///////////////////////////////////////////////////////////////////////////////
template <typename T>
class internal_deque: public linear_memory_base<internal_deque<T> > {
	array<T> m_elements;
	memory_size_type m_first, m_size;
public:
	///////////////////////////////////////////////////////////////////////////
	/// \copybrief linear_memory_structure_doc::memory_coefficient()
	/// \copydetails linear_memory_structure_doc::memory_coefficient()
	///////////////////////////////////////////////////////////////////////////
	static double memory_coefficient() {
		return array<T>::memory_coefficient();
	}

	///////////////////////////////////////////////////////////////////////////
	/// \copybrief linear_memory_structure_doc::memory_overhead()
	/// \copydetails linear_memory_structure_doc::memory_overhead()
	///////////////////////////////////////////////////////////////////////////
	static double memory_overhead() {
		return array<T>::memory_overhead() - sizeof(array<T>) + sizeof(internal_deque);
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Construct a deque.
	///
	/// \param size The number of pushes supported between calls to clear and
	/// resize.
	///////////////////////////////////////////////////////////////////////////
	internal_deque(size_t size = 0): m_first(0), m_size(0) {
		m_elements.resize(size);
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Resize the deque; all data is lost.
	///
	/// \param size The number of elements to contain.
	///////////////////////////////////////////////////////////////////////////
	void resize(size_t size = 0) {
		m_elements.resize(size);
		m_first = m_size = 0;
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Return the first item in the deque.
	///////////////////////////////////////////////////////////////////////////
	T front() {
		tp_assert(!empty(), "front() on an empty deque");
		return m_elements[m_first];
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Return the last item in the deque.
	///////////////////////////////////////////////////////////////////////////
	T back() {
		tp_assert(!empty(), "back() on an empty deque");
		memory_size_type index = (m_first + m_size) % m_elements.size();
		if(index == 0) // avoid overflow
			index = m_elements.size() - 1;
		else
			--index;

		return m_elements[index];
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Add an element to the back of the deque.
	///
	/// \param val The element to add.
	///////////////////////////////////////////////////////////////////////////
	void push_back(T val) {
		tp_assert(!full(), "push_back() on a full deque");
		m_elements[(m_first + m_size) % m_elements.size()] = val;
		++m_size;
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Remove an element from the back of the deque.
	///////////////////////////////////////////////////////////////////////////
	void pop_back(){
		tp_assert(!empty(), "pop_back() on an empty deque");
		--m_size;
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Add an element to the front of the deque.
	///
	/// \param val The element to add.
	///////////////////////////////////////////////////////////////////////////
	void push_front(T val) {
		tp_assert(!full(), "push_front() on a full deque");
		if(m_first == 0) // avoid overflow
			m_first = m_elements.size() -1;
		else
			m_first = (m_first-1) % m_elements.size();
		++m_size;
		m_elements[m_first] = val;
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Remove an element from the front of the deque.
	///////////////////////////////////////////////////////////////////////////
	void pop_front(){
		tp_assert(!empty(), "pop_front() on an empty deque");
		m_first = (m_first+1) % m_elements.size();
		--m_size;
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Check if the deque is empty.
	/// \return true if the deque is empty, otherwise false.
	///////////////////////////////////////////////////////////////////////////
	bool empty() const {
		return size() == 0;
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Check if the deque is full.
	/// \return true if the deque is full, otherwise false.
	///////////////////////////////////////////////////////////////////////////
	bool full() const {
		return size() == m_elements.size();
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Return the number of elements in the deque.
	/// \return The number of elements in the deque.
	///////////////////////////////////////////////////////////////////////////
	size_t size() const {
		return m_size;
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Clear the deque of all elements.
	/// After this call, the deque again supports the number of pushes passed
	/// to the constructor or resize.
	///////////////////////////////////////////////////////////////////////////
	void clear() {
		m_first = m_size = 0;
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Element access.
	///////////////////////////////////////////////////////////////////////////////
	T & operator[](memory_size_type i) {
		return m_elements[(m_first+i) % m_elements.size()];
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Element access.
	///////////////////////////////////////////////////////////////////////////////
	const T & operator[](memory_size_type i) const {
		return m_elements[(m_first+i) % m_elements.size()];
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Element access with bounds checking
	///////////////////////////////////////////////////////////////////////////////
	T & at(memory_size_type i) {
		tp_assert(i < size(), "index out of bounds");
		return m_elements[(m_first+i) % m_elements.size()];
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Element access with bounds checking
	///////////////////////////////////////////////////////////////////////////////
	const T & at(memory_size_type i) const {
		tp_assert(i < size(), "index out of bounds");
		return m_elements[(m_first+i) % m_elements.size()];
	}
}; // class internal_deque

} // namespace tpie
#endif //__TPIE_INTERNAL_DEQUE_H__

