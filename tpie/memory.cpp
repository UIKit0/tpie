// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :
//
// Copyright 2011, The TPIE development team
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

#include "memory.h"
#include <iostream>
#include <sstream>
#include "tpie_log.h"
#include "static_string_stream.h"

namespace tpie {

memory_manager * mm = 0;

memory_manager::memory_manager(): m_used(0), m_limit(0), m_maxExceeded(0), m_enforce(ENFORCE_WARN) {}

size_t memory_manager::available() const throw() {
	size_t used = m_used;
	size_t limit = m_limit;
	if (used < limit) return used-limit;
	return 0;
}

void memory_manager::register_allocation(size_t bytes) {
	boost::mutex::scoped_lock lock(m_mutex);

	switch(m_enforce) {
	case ENFORCE_IGNORE:
		m_used += bytes;
		break;
	case ENFORCE_THROW:
		if (m_used + bytes > m_limit) {
			lock.unlock();
			std::stringstream ss;
			ss << "Memory allocation error, memory limit exceeded. "
			   <<"Allocation request \""
			   << bytes
			   << "\" plus previous allocation \""
			   << m_used
			   << "\" exceeds user-defined limit \""
			   << m_limit
			   << "\"";
			throw out_of_memory_error(ss.str().c_str());
		}
		m_used += bytes;
		break;
	case ENFORCE_WARN:
		m_used += bytes;
		if (m_used > m_limit && m_used - m_limit > m_maxExceeded) {
			m_maxExceeded = m_used - m_limit;
			lock.unlock();
			log_warning() << "Memory limit exceeded by " << m_maxExceeded 
						  << " bytes, while trying to allocate " << bytes << " bytes" << std::endl;
		}
	};
}

void memory_manager::register_deallocation(size_t bytes) {
	boost::mutex::scoped_lock lock(m_mutex);
#ifndef TPIE_NDEBUG
	if (bytes > m_used)
		log_error() << "Error in deallocation, trying to deallocate " << bytes << " bytes, while only " <<
			m_used << " where allocated" << std::endl;
#endif
	m_used -= bytes;
}


void memory_manager::set_limit(size_t new_limit) {
	boost::mutex::scoped_lock lock(m_mutex);
	m_limit = new_limit;
}

void memory_manager::set_enforcement(enforce_t e) {
	boost::mutex::scoped_lock lock(m_mutex);
	m_enforce = e;
}

struct log_flusher {
	tpie::static_string_stream buf;
	~log_flusher() {
		if(strlen(buf.c_str())) {
			tpie::log_debug() << buf.c_str();
			tpie::log_debug().flush();
		}
	}
};

std::pair<uint8_t *, size_t> memory_manager::__allocate_consecutive(size_t upper_bound, size_t granularity) {
	log_flusher lf;

	boost::mutex::scoped_lock lock(m_mutex);

	size_t high=available()/granularity;
	if (upper_bound != 0) 
		high=std::min(upper_bound/granularity, high);
	size_t low=1;
	uint8_t * res;

	//first check quickly if we can get "high" bytes of memory
	//directly.
	try {
		res = new uint8_t[high*granularity];
		m_used += high*granularity;
		return std::make_pair(res, high*granularity);
	} catch (std::bad_alloc) {
		lf.buf << "Failed to get " << (high*granularity)/(1024*1024) << " megabytes of memory. "
		 	   << "Performing binary search to find largest amount "
			   << "of memory available. This might take a few moments.\n";
	}

	//perform a binary search in [low,high] for highest possible 
	//memory allocation
	size_t best=0;
	do {
		size_t mid = static_cast<size_t>((static_cast<uint64_t>(low)+high)/2);

		lf.buf << "Search area is  [" << low*granularity << "," << high*granularity << "]"
			   << " query amount is: " << mid*granularity << ":\n";

		//try to allocate "mid" bytes of memory
		//std::new throws an exception if memory allocation fails
		try {
			uint8_t* mem = new uint8_t[mid*granularity];
			low = mid+1;
			best=mid*granularity;
			delete[] mem;
		} catch (std::bad_alloc) {
			high = mid-1;
			lf.buf << "   failed.\n";
		}
	} while (high >= low);
	lf.buf << "- - - - - - - END MEMORY SEARCH - - - - - -\n";	

	res = new uint8_t[best];
	m_used += best;
	return std::make_pair(res, best);
}

void init_memory_manager() {
	mm = new memory_manager();
}

void finish_memory_manager() {
	delete mm;
	mm = 0;
}

memory_manager & get_memory_manager() {return * mm;}


size_t consecutive_memory_available(size_t granularity) {
	std::pair<uint8_t *, size_t> r = get_memory_manager().__allocate_consecutive(0, granularity);
	tpie_delete_array(r.first, r.second);
	return r.second;
}


} //namespace tpieg