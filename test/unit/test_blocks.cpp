// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; c-file-style: "stroustrup"; -*-
// vi:set ts=4 sts=4 sw=4 noet cino+=(0 :
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

#include "common.h"
#include <tpie/blocks/b_tree.h>
#include <tpie/prime.h>
#include <tpie/pipelining/b_tree.h>
#include <numeric>
#include <boost/random/linear_congruential.hpp>
#include <set>
#include <boost/function_output_iterator.hpp>

typedef size_t key_type;
typedef tpie::blocks::b_tree<key_type> tree_type;
typedef tpie::blocks::b_tree_builder<key_type> builder_type;

bool b_tree_test() {
	bool result = true;
	tree_type t;
	t.open();
	for (key_type i = 0; i < 100; ++i) {
		t.insert((3*i) % 100);
		if (!t.count((i/2*3) % 100)) {
			tpie::log_error() << "Missing element " << (i/2*3)%100 << " from B tree" << std::endl;
			result = false;
		}
	}
	return result;
}

bool b_tree_test_2(key_type items) {
	key_type p = static_cast<key_type>(tpie::next_prime(static_cast<size_t>(items+1)));
	tpie::log_debug() << "Generating items " << p << "*i%" << items << " for i in [0," << items << ")" << std::endl;
	tree_type t;
	t.open();
	for (key_type i = 0; i < items; ++i) {
		try {
			t.insert(p*i%items);
		} catch (...) {
			tpie::log_error() << "Exception after " << i << " insertions" << std::endl;
			throw;
		}
	}
	std::vector<key_type> output;
	output.reserve(items);
	t.in_order_dump(std::back_inserter(output));

	if (output.size() != items) {
		tpie::log_error() << "B tree dump output incorrect no. of items" << std::endl;
		return false;
	}

	for (key_type i = 0; i < items; ++i) {
		if (output[i] != i) {
			tpie::log_error() << "B tree dump incorrect @ " << i << std::endl;
			return false;
		}
	}
	return true;
}

struct average_augment {
	int sum;
	size_t count;

	float average() const {
		return static_cast<float>(sum)/count;
	}

	bool operator==(const average_augment & other) {
		return sum == other.sum && count == other.count;
	}
};

std::ostream & operator<<(std::ostream & out, const average_augment & a) {
	return out << "[sum=" << a.sum << "; count=" << a.count << "]";
}

struct average_augmentor {
	average_augment operator()(int * first, int * last) const {
		average_augment res;
		res.sum = std::accumulate(first, last, 0);
		res.count = static_cast<size_t>(last-first);

		return res;
	}

	average_augment operator()(average_augment * first, average_augment * last) const {
		average_augment res;
		res.sum = 0;
		res.count = 0;

		for(average_augment * i = first; i != last; ++i) {
			res.sum += i->sum;
			res.count += i->count;
		}

		return res;
	}
};

bool b_tree_augmented_test() {
	typedef tpie::blocks::b_tree<int, int, std::less<int>, tpie::blocks::identity_key_extract<int>, average_augment , average_augmentor> tree_t;

	size_t items = 10000;
	size_t mod = 100;
	boost::rand48 random(42);

	tree_t tree1;
	std::multiset<int> tree2;

	tree1.open();

	for(key_type i = 0; i < items; ++i) {
		try {
			size_t n = random() % mod;
			tree1.insert(n);
			tree2.insert(n);

			size_t a = random() % mod;
			size_t b = random() % mod; b += a;
			// calculate the correct average
			average_augment correct;
			correct.sum = 0;
			correct.count = 0;

			std::set<int>::iterator begin = tree2.lower_bound(a);
			std::set<int>::iterator end = tree2.upper_bound(b);

			for(std::set<int>::iterator i = begin; i != end; ++i) {
				correct.sum += *i;
				++correct.count;
			}

			average_augment res = tree1.augment(a, b);

			TEST_ENSURE_EQUALITY(correct, res, "the b_tree did not output the correct result.")
		} catch (...) {
			std::stringstream ss;
			ss << "Exception after " << i << " insertions";
			TEST_FAIL(ss.str());
		}
	}

	return true;
}

class verify_iterator {
public:
	verify_iterator(size_t * outputSize, bool * error, key_type * i, key_type a, key_type step, key_type b)
		: outputSize(outputSize)
		, error(error)
		, i(i)
		, a(a)
		, step(step)
		, b(b)
	{
	}

	void operator++() const {}

	const verify_iterator & operator*() const {
		return *this;
	}

	const verify_iterator & operator=(key_type v) const {
		++*outputSize;
		if (*error) return *this;
		if (v != *i) {
			tpie::log_error() << "B tree dump incorrect: got " << v << ", expected " << *i << std::endl;
			*error = true;
		} else {
			*i += step;
		}
		return *this;
	}

private:
	size_t * outputSize;
	bool * error;
	key_type * i;
	key_type a;
	key_type step;
	key_type b;
};

bool verify_tree(tree_type & t, key_type a, key_type step, key_type b) {
	size_t items = (a == b) ? 0 : ((b-a+step-1)/step);
	size_t outputSize = 0;
	key_type i = a;
	bool error = false;
	t.in_order_dump(verify_iterator(&outputSize, &error, &i, a, step, b));
	if (outputSize != items) {
		tpie::log_error()
			<< "B tree dump output incorrect no. of items\n"
			<< "Expected " << items << "; got " << outputSize << '\n'
			<< std::flush;
		return false;
	}
	return !error;
}

bool b_tree_erase_test(key_type items, size_t fanout) {
	tree_type t;
	t.open();
	if (fanout != 0) {
		tpie::blocks::b_tree_parameters params = t.get_parameters();
		params.nodeMax = params.leafMax = fanout;
		params.nodeMin = params.leafMin = (fanout + 3) / 4;
		t.set_parameters(params);
	}
	{
		builder_type b(t);
		for (key_type i = 0; i < items; ++i) {
			b.push(i);
		}
		b.end();
	}
	if (!verify_tree(t, 0, 1, items)) return false;
	for (key_type i = 0; i < items; i += 2) {
		t.erase(i);
	}
	if (!verify_tree(t, 1, 2, items)) return false;
	for (key_type i = 0; i < items; i += 2) {
		t.insert(i);
	}
	if (!verify_tree(t, 0, 1, items)) return false;
	for (key_type i = 0; i < items; ++i) {
		t.erase(i);
	}
	if (!verify_tree(t, 0, 0, 0)) return false;
	return true;
}

bool b_tree_builder_test(key_type n) {
	tree_type t;
	t.open();
	{
		builder_type builder(t);
		for (key_type i = 0; i < n; ++i) builder.push(i);
		builder.end();
	}
	if (!verify_tree(t, 0, 1, n)) return false;
	return true;
}

template <typename dest_t>
class iota_t : public tpie::pipelining::node {
	key_type n;
	dest_t dest;
public:
	iota_t(dest_t dest, key_type n)
		: n(n)
		, dest(dest)
	{
		this->set_steps(n);
		this->add_push_destination(dest);
		this->set_name("Iota");
	}

	virtual void go() {
		for (key_type i = 0; i < n; ++i) {
			dest.push(i);
			this->step();
		}
	}
};

tpie::pipelining::pipe_begin<tpie::pipelining::factory_1<iota_t, key_type> >
iota(key_type n) {
	return tpie::pipelining::factory_1<iota_t, key_type>(n);
}

bool b_tree_builder_pipe_test(key_type n) {
	tree_type t;
	t.open();
	{
		builder_type builder(t);
		tpie::pipelining::pipeline p =
			iota(n) | tpie::pipelining::b_tree_builder(builder);
		p.plot(tpie::log_debug());
		p();
	}
	if (!verify_tree(t, 0, 1, n)) return false;
	return true;
}

struct range_report_tester {
	range_report_tester(size_t a, size_t b, bool & correct, size_t & count) : a(a), b(b), correct(correct), count(count) {}

	void operator()(const size_t & x) {
		++count;
		tpie::log_debug() << "Range reporter given element " << x << std::endl;
		if(!(a <= x && x <= b)) {
			correct = false;
			tpie::log_error() << "Range reporter given element" << x << " which is outside the range [" << a << "; " << b << "]." << std::endl;
		}
	}
private:
	size_t a, b;
	bool & correct;
	size_t & count;
};

bool b_tree_range_report_test() {
	size_t items = 10000;
	boost::rand48 random(42);

	tpie::blocks::b_tree<size_t> tree1;
	std::multiset<size_t> tree2;

	tree1.open();

	for(key_type i = 0; i < items; ++i) {
		try {
			size_t n = random() % items;
			size_t a = random() % items;
			size_t b = random() % items;
			if(a > b) std::swap(a, b);

			tree1.insert(n);
			tree2.insert(n);

			bool correct = true;
			size_t count = 0;
			size_t correctCount = std::distance(tree2.lower_bound(a), tree2.upper_bound(b));

			range_report_tester tester(a, b, correct, count);
			tree1.range_report(a, b, boost::make_function_output_iterator(tester));

			TEST_ENSURE(correct, "range_report reported an element outside the range");

			if(count != correctCount) {
				tpie::log_error() << "range_report did not report the correct number of elements" << std::endl;
				tpie::log_error() << "Element: " << i << std::endl;
				tpie::log_error() << "Range: [" << a << "; " << b << "]" << std::endl;
				tpie::log_error() << count << " != " << correctCount << std::endl;
			}
		} catch (...) {
			std::stringstream ss;
			ss << "Exception after " << i << " insertions";
			TEST_FAIL(ss.str());
		}
	}


	return true;
}

int main(int argc, char ** argv) {
	return tpie::tests(argc, argv)
	.test(b_tree_test, "b_tree")
	.test(b_tree_test_2, "b_tree_2", "n", static_cast<key_type>(1000))
	.test(b_tree_augmented_test, "b_tree_augmented")
	.test(b_tree_erase_test, "b_tree_erase",
		  "n", static_cast<key_type>(1000),
		  "fanout", static_cast<size_t>(0))
	.test(b_tree_builder_test, "b_tree_builder", "n", static_cast<key_type>(1000))
	.test(b_tree_builder_pipe_test, "b_tree_builder_pipe", "n", static_cast<key_type>(1000))
	.test(b_tree_range_report_test, "b_tree_range_report")
	;
}
