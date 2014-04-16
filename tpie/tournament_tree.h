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
#ifndef __TPIE_TOURNAMENT_TREE_H__
#define __TPIE_TOURNAMENT_TREE_H__

///////////////////////////////////////////////////////////////////////////////
/// \file internal_queue.h
/// Generic internal queue with known memory requirements.
///////////////////////////////////////////////////////////////////////////////
#include <tpie/array.h>
#include <tpie/util.h>
#include <tpie/tpie_assert.h>

namespace tpie {

namespace bits {

template<typename T, typename pred_t>
class tournament_tree {
public:
	tournament_tree(T begin, T end, pred_t pred) : m_elements(end-begin-1), m_begin(begin), m_end(end), m_pred(pred) {
		const memory_size_type n = end-begin;
		const memory_size_type k = n/2-1;
		for(memory_size_type i = 0; k+i < n-1; ++i) {
			if(!m_pred(*(m_begin+2*i+1), *(m_begin+2*i))) {
				m_elements[i+k] = 2*i;
			}
			else {
				m_elements[i+k] = 2*i+1;
			}
		}

		for(memory_size_type i = k; i > 0;) {
			--i;
			if(!m_pred(*(m_begin+m_elements[2*i+2]), *(m_begin+m_elements[2*i+1]))) {
				m_elements[i] = m_elements[2*i+1];
			}
			else {
				m_elements[i] = m_elements[2*i+2];
			}
		}
	}

	tournament_tree(pred_t pred) : m_pred(pred) {}

	memory_size_type top() {
		return m_elements[0];
	}

	void update_key(memory_size_type i) {
		const memory_size_type n = m_end-m_begin;
		const memory_size_type parent = i/2 + n/2 - 1;
		if(!m_pred(*(m_begin+(i^1)), *(m_begin + i))) {
			m_elements[parent] = i;
		}
		else {
			m_elements[parent] = i^1;
		}

		memory_size_type j = parent;
		while(j != 0) {
			j = (j-1) / 2;
			if(!m_pred(*(m_begin+m_elements[2*j+2]), *(m_begin+m_elements[2*j+1]))) {
				m_elements[j] = m_elements[2*j+1];
			}
			else {
				m_elements[j] = m_elements[2*j+2];
			}
		}
	}
private:
	tpie::array<memory_size_type> m_elements;
	T m_begin;
	T m_end;
	pred_t m_pred;
};

} // namespace bits

} // namespace tpie
#endif //__TPIE_TOURNAMENT_TREE_H__

