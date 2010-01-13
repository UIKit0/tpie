 // -*- mode: c++; tab-width: 4; indent-tabs-mode: t; c-file-style: "stroustrup"; -*-
// vi:set ts=4 sts=4 sw=4 noet :
// Copyright 2009, The TPIE development team

// This file is part of TPIE.

// TPIE is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.

// TPIE is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with TPIE.  If not, see <http://www.gnu.org/licenses/>
#include "../app_config.h"
#include <tpie/portability.h>
#include <cstring>
// #include <tpie/stream/posix_bte.h>
// #include <tpie/stream/concepts.h>
// #include <tpie/stream/exception.h>
// #include <tpie/stream/fd_file.h>
#include <tpie/file_stream.h>
#include <tpie/file.h>
#include <tpie/util.h>
#include <tpie/file_accessor/stdio.h>

using namespace std;
using namespace tpie;

#define ERR(x) {cerr << x << endl; exit(1);}

template <typename T>
void test_file_accessor() {
 	remove("tmp");
 	{
 		int d=42;

 		T x;
 		x.open("tmp", false, true, sizeof(int), sizeof(int));
 		if (x.size() != 0) ERR("New stream has wrong size");
 		if (x.path() != "tmp") ERR("Wrong path");
		
 		x.write(&d, 0, 1);
 		x.write(&d, 1, 1);
		
		int ud=314;
		x.write_user_data(&ud);

 		try {
 			x.read(&d, 0, 1);
 			ERR("Read should faild");
 		} catch(io_exception &) {
 			//Do nothing
 		}
		
 		if (x.size() != 2) ERR("Wrong size");
 		x.close();
 	}

 	try {
 		T x;
 		x.open("tmp", true, false, sizeof(int)+1, sizeof(int));
 		ERR("Opened file with wrong item size");
 	} catch(invalid_file_exception&) {
 		//Do nothing
 	}

 	try {
 		T x;
 		x.open("tmp", true, false, sizeof(int), 0);
 		ERR("Opened file with wrong user data size");
 	} catch(invalid_file_exception&) {
 		//Do nothing
 	}

 	{
 		int d;
 		T x;
 		x.open("tmp", true, true, sizeof(int), sizeof(int));
 		if (x.read(&d, 1, 1) != 1 || d != 42) ERR("Read failed");
 		d=12;
 		x.write(&d, 1, 1);
 		x.write(&d, 2, 1);
 		if (x.read(&d, 0, 1) != 1 || d != 42) ERR("Read failed");
 		if (x.read(&d, 1, 1) != 1 || d != 12) ERR("Read failed");
 		if (x.read(&d, 2, 1) != 1 || d != 12) ERR("Read failed");
		int ud;
		x.read_user_data(&ud);

		if (ud != 314) ERR("Wrong user data");
 		if (x.size() != 3) ERR("Wrong size");
 		x.close();
 	}
	
 	{
 		T x;
 		x.open("tmp", true, false, sizeof(int), sizeof(int) );
 		try {
 			int d=44;
 			x.write(&d, 0, 1);
 			ERR("Write should faild");
 		} catch(io_exception &) {
 			//Do nothing
 		}
 		int d;
 		if (x.read(&d, 0, 1) != 1 || d != 42) ERR("Read failed");
 		if (x.read(&d, 1, 1) != 1 || d != 12) ERR("Read failed");
 		if (x.read(&d, 2, 1) != 1 || d != 12) ERR("Read failed");
 		x.close();
 	}
 	remove("tmp");
}

int main(int argc, char ** argv) {
 	if (argc == 2 && !strcmp(argv[1], "file_accessor_stdio")) {
 		test_file_accessor<file_accessor::stdio>();
 	} else if (argc == 2 && !strcmp(argv[1], "file_stream")) {
		///First a simple test
		remove("tmp");
		{
			file_stream<int> stream;
			stream.open("tmp", file_base::write, sizeof(int));

			stream.write_user_data<int>(42);
			if (stream.size() != 0) ERR("size failed(1)");
			for(int i=0; i < 40; ++i)
				stream.write((i*8209)%8273);
 			if (stream.size() != 40) ERR("size failed(2)");
 			stream.close();
		}

		{
			file_stream<int> stream;
			stream.open("tmp", file_base::read, sizeof(int));
			if (stream.size() != 40) ERR("size failed(3)");
			for(int i=0; i< 40; ++i) {
				if (stream.has_more() == false) ERR("has_more failed");
				if (stream.read() != (i*8209)%8273) ERR("read failed");
			}
			if (stream.has_more() == true) ERR("has_more failed (2)");
			try {
				int r =stream.read();
				unused(r);
				ERR("read did not fail as expected");
			} catch(end_of_stream_exception &) {
				//Do nothing
			}		
			
			int y;
			stream.read_user_data<int>(y);
			if (y != 42) ERR("read did not fail as expected");

			stream.close();
		}
	} else if (argc == 2 && !strcmp(argv[1], "substreams")) {
		tpie::remove("tmp");
		file<int> f;
		f.open("tmp", file_base::read_write);
 			
		tpie::seed_random(1234);
		const int cnt=4;
		const int size=128;
		typedef file<int>::stream stream_t;
		stream_t ** streams = new stream_t*[cnt];
		for(int i=0; i < cnt; ++i) 
			streams[i] = new stream_t(f);
		
		int content[size];
		for(int i=0; i < size; ++i) {
			content[i] = tpie::random();
			streams[0]->write(content[i]);
		}
		
		for(int i=0; i < 200000; ++i ) {
			int l=tpie::random()%size;
			int s=tpie::random()%cnt;
			streams[s]->seek(l);
			if (tpie::random() % 2 == 0) {
				if (streams[s]->read() != content[l]) ERR("read failed(2)");
			} else {
				content[l] = tpie::random();
				streams[s]->write(content[l]);
			}
		}
		for(int i=0; i < cnt; ++i)
			delete streams[i];
		delete[] streams;
	} else {
		return 1;
	}
	return 0;
}