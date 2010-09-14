// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :
// Copyright 2010, The TPIE development team
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

#ifndef __TPIE_CONCEPTS_H__
#define __TPIE_CONCEPTS_H__
#include <tpie/types.h>
#include <boost/concept_check.hpp>
namespace tpie {

/////////////////////////////////////////////////////////
/// \file concepts.h
/// \brief Implementaton boost concept checkers
/////////////////////////////////////////////////////////

namespace concepts {

template <class T>
class memory_calculatable {
public:
	BOOST_CONCEPT_USAGE(memory_calculatable) {
		memory_size_type n = T::memory(42);
		unused(n);
	}
};

}

/////////////////////////////////////////////////////////
/// \brief Check if a structure adhears to the linear_memory_structure concept
/// 
/// \sa linear_memory_structure_doc
/////////////////////////////////////////////////////////
template <typename T>
struct linear_memory_structure_concept {
  BOOST_CONCEPT_USAGE(linear_memory_structure_concept) {
	  (double)T::memory_coefficient();
	  (double)T::memory_overhead();
  }
};

}
#endif //__TPIE_CONCEPTS_H__

