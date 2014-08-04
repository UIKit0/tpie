// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; c-file-style: "stroustrup"; -*-
// vi:set ts=4 sts=4 sw=4 noet cino=(0 :
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

// block_collection usage test

#include "common.h"
#include <tpie/tpie.h>
#include <tpie/blocks/block_collection.h>
#include <tpie/tempname.h>
#include <vector>
#include <tpie/file_accessor/file_accessor.h>

using namespace tpie;
using namespace tpie::blocks;

const stream_size_type max_block_size = 5 * 1024 * 1024;

memory_size_type random(memory_size_type i) {
	return 179424673 * i + 15485863;
}

bool basic() {
	std::vector<block_handle> blocks;

	temp_file file;
	block_collection collection(file.path(), true);

	// write 20 twenty blocks of random sizes
	for(char i = 0; i < 20; ++i) {
		stream_size_type size = random(i) % max_block_size + 1;

		block_handle handle = collection.get_free_block(size);

		TEST_ENSURE(handle.size >= size, "The returned block size is too small.");

		block b(handle.size);

		for(block::iterator j = b.begin(); j != b.end(); ++j)
			*j = i;

		collection.write_block(handle, b);
		blocks.push_back(handle);
	}

	// close and reopen the collection
	collection.close();
	collection.open(file.path(), true);

	// verify the content of the 20 blocks
	for(char i = 0; i < 20; ++i) {
		block_handle handle = blocks[i];

		block b;
		collection.read_block(handle, b);

		TEST_ENSURE_EQUALITY(handle.size, b.size(), "The block size should be equal to the handle size");

		for(block::iterator j = b.begin(); j != b.end(); ++j)
			TEST_ENSURE_EQUALITY((int) *j, (int) i, "the content of the returned block is not correct");
	}

	collection.close();

	return true;
}

int main(int argc, char **argv) {
	return tpie::tests(argc, argv)
		.test(basic, "basic");

}