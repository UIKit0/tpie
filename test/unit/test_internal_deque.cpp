// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; c-file-style: "stroustrup"; -*-
// vi:set ts=4 sts=4 sw=4 noet :
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

#include "common.h"

#include <tpie/internal_deque.h>

using namespace tpie;

size_t random(size_t i) {
	return (i * 104729) % 2251;
}

bool basic_test() {
	internal_deque<size_t> q(52);

	TEST_ENSURE(!q.full(), "full() should be false");
	TEST_ENSURE(q.empty(), "empty() should be true");

	for(size_t i=0; i < 52; ++i) {
		size_t j = random(i);
		if(i % 2 == 0) {
			q.push_back(j);
			TEST_ENSURE_EQUALITY(q.back(), j, "back() did not return the correct value")
		}
		else {
			q.push_front((i * 104729) % 2251);
			TEST_ENSURE_EQUALITY(q.front(), j, "front() did not return the correct value")
		}

		TEST_ENSURE_EQUALITY(i, q.size()-1, "size() did not return the correct value");
		if(i == 0)
			TEST_ENSURE_EQUALITY(q.front(), q.back(), "front() and back() should be equal when size == 1")
	}

	TEST_ENSURE(q.full(), "full() should be true");
	TEST_ENSURE(!q.empty(), "empty() should be false");

	for(size_t i = 0; i < 52; ++i) {
		TEST_ENSURE_EQUALITY(q.size(), 52-i, "size() did not return the correct result");

		size_t k = 52-i-1;
		if(k % 2 == 0) {
			TEST_ENSURE_EQUALITY(q.back(), random(k), "back() did not return the correct result");
			q.pop_back();
		}
		else {
			TEST_ENSURE_EQUALITY(q.front(), random(k), "front() did not return the correct result");
			q.pop_front();
		}
	}

	TEST_ENSURE(q.empty(), "empty() did not return the correct result");
	return true;
}

class deque_memory_test: public memory_test {
public:
	internal_deque<size_t> * a;
	virtual void alloc() {a = tpie_new<internal_deque<size_t> >(123456);}
	virtual void free() {tpie_delete(a);}
	virtual size_type claimed_size() {return static_cast<size_type>(internal_deque<size_t>::memory_usage(123456));}
};

int main(int argc, char **argv) {
	return tpie::tests(argc, argv)
		.test(basic_test, "basic")
		.test(deque_memory_test(), "memory");
}


