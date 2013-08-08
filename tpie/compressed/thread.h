// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; c-file-style: "stroustrup"; -*-
// vi:set ts=4 sts=4 sw=4 noet :
// Copyright 2013, The TPIE development team
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

#ifndef TPIE_COMPRESSED_THREAD_H
#define TPIE_COMPRESSED_THREAD_H

#include <boost/thread.hpp>
#include <tpie/array.h>
#include <tpie/file_accessor/file_accessor.h>
#include <tpie/file_accessor/byte_stream_accessor.h>

namespace tpie {

class compressor_thread;

void init_compressor();
void finish_compressor();

compressor_thread & the_compressor_thread();

class compressor_buffer {
private:
	typedef array<char> storage_t;

	storage_t m_storage;
	memory_size_type m_size;

public:
	compressor_buffer(memory_size_type capacity)
		: m_storage(capacity)
		, m_size(0)
	{
	}

	~compressor_buffer() {
	}

	char * get() {
		return m_storage.get();
	}

	const char * get() const {
		return m_storage.get();
	}

	memory_size_type size() const {
		return m_size;
	}

	memory_size_type capacity() const {
		return m_storage.size();
	}

	void set_size(memory_size_type size) {
		m_size = size;
	}

	void set_capacity(memory_size_type capacity) {
		m_storage.resize(capacity);
		m_size = 0;
	}
};

class compressor_thread_lock;

#ifdef __GNUC__
class __attribute__((__may_alias__)) read_request;
class __attribute__((__may_alias__)) write_request;
#endif // __GNUC__

class read_request {
	struct state {
		bool done;
		bool endOfStream;
		stream_size_type nextReadOffset;
		stream_size_type nextBlockSize;
	};
public:
	typedef boost::shared_ptr<compressor_buffer> buffer_t;
	typedef file_accessor::byte_stream_accessor<default_raw_file_accessor> file_accessor_t;
	typedef boost::condition_variable condition_t;

	read_request(buffer_t buffer,
				 file_accessor_t * fileAccessor,
				 stream_size_type readOffset,
				 stream_size_type blockSize,
				 boost::condition_variable & cond)
		: m_buffer(buffer)
		, m_fileAccessor(fileAccessor)
		, m_readOffset(readOffset)
		, m_blockSize(blockSize)
		, m_cond(cond)
		, m_state(new state)
	{
		m_state->done = false;
		m_state->endOfStream = false;
		m_state->nextReadOffset = 0;
	}

	inline void wait(compressor_thread_lock & lock);

	void notify() {
		m_cond.notify_one();
	}

	buffer_t buffer() {
		return m_buffer;
	}

	bool done() const {
		return m_state->done;
	}

	void set_done() {
		m_state->done = true;
	}

	bool end_of_stream() const {
		return m_state->endOfStream;
	}

	void set_end_of_stream() {
		m_state->endOfStream = true;
	}

	file_accessor_t & file_accessor() {
		return *m_fileAccessor;
	}

	stream_size_type read_offset() {
		return m_readOffset;
	}

	stream_size_type block_size() {
		return m_blockSize;
	}

	stream_size_type next_read_offset() {
		return m_state->nextReadOffset;
	}

	void set_next_read_offset(stream_size_type o) {
		m_state->nextReadOffset = o;
	}

	stream_size_type next_block_size() {
		return m_state->nextBlockSize;
	}

	void set_next_block_size(stream_size_type o) {
		m_state->nextBlockSize = o;
	}

private:
	buffer_t m_buffer;
	file_accessor_t * m_fileAccessor;
	/** If readOffset is zero, the next block to read is the first block and its size is not known.
	 * In that case, the size of the first block is the first eight bytes, and the first block begins
	 * after those eight bytes.
	 * If readOffset and blockSize are both non-zero, the next block begins at the given offset
	 * and has the given size.
	 * Otherwise, if readOffset is non-zero and blockSize is zero, we have reached the end of
	 * the stream.
	 */
	const stream_size_type m_readOffset;
	const stream_size_type m_blockSize;
	condition_t & m_cond;
	boost::shared_ptr<state> m_state;
};

class write_request {
public:
	typedef boost::shared_ptr<compressor_buffer> buffer_t;
	typedef file_accessor::byte_stream_accessor<default_raw_file_accessor> file_accessor_t;

	write_request(const buffer_t & buffer, file_accessor_t * fileAccessor, memory_size_type blockItems)
		: m_buffer(buffer)
		, m_fileAccessor(fileAccessor)
		, m_blockItems(blockItems)
	{
	}

	file_accessor_t & file_accessor() {
		return *m_fileAccessor;
	}

	buffer_t buffer() {
		return m_buffer;
	}

	memory_size_type block_items() {
		return m_blockItems;
	}

private:
	buffer_t m_buffer;
	file_accessor_t * m_fileAccessor;
	const memory_size_type m_blockItems;
};

class compressor_request_kind {
public:
	enum type {
		NONE,
		READ,
		WRITE
	};

private:
	compressor_request_kind() /*= delete*/;
	compressor_request_kind(const compressor_request_kind &) /*= delete*/;
	~compressor_request_kind() /*= delete*/;
};

///////////////////////////////////////////////////////////////////////////////
/// \brief  Tagged union containing either a read_request or a write_request.
///
/// In C++11, this can be implemented more elegantly since the C++11 standard
/// allows unions with member types having constructors and destructors.
/// In C++03, we have to emulate this with a char buffer that is large enough.
///
/// To turn the tagged union into a T, call set_T(ctor params) which returns a
/// reference to T. When the tagged union is a T (check with kind()), use
/// get_T() to get the reference to T.
///////////////////////////////////////////////////////////////////////////////
class compressor_request {
public:
	compressor_request()
		: m_kind(compressor_request_kind::NONE)
	{
	}

	~compressor_request() {
		destruct();
	}

	compressor_request(const compressor_request & other)
		: m_kind(compressor_request_kind::NONE)
	{
		switch (other.kind()) {
			case compressor_request_kind::NONE:
				break;
			case compressor_request_kind::READ:
				set_read_request(other.get_read_request());
				break;
			case compressor_request_kind::WRITE:
				set_write_request(other.get_write_request());
				break;
		}
	}

	read_request & set_read_request(const read_request::buffer_t & buffer,
									read_request::file_accessor_t * fileAccessor,
									stream_size_type readOffset,
									stream_size_type blockSize,
									boost::condition_variable & cond)
	{
		destruct();
		m_kind = compressor_request_kind::READ;
		return *new (m_payload) read_request(buffer, fileAccessor, readOffset, blockSize, cond);
	}

	read_request & set_read_request(const read_request & other) {
		destruct();
		m_kind = compressor_request_kind::READ;
		return *new (m_payload) read_request(other);
	}

	write_request & set_write_request(const write_request::buffer_t & buffer,
									  write_request::file_accessor_t * fileAccessor,
									  memory_size_type blockItems)
	{
		destruct();
		m_kind = compressor_request_kind::WRITE;
		return *new (m_payload) write_request(buffer, fileAccessor, blockItems);
	}

	write_request & set_write_request(const write_request & other) {
		destruct();
		m_kind = compressor_request_kind::WRITE;
		return *new (m_payload) write_request(other);
	}

	// Precondition: kind() == READ
	read_request & get_read_request() {
		return *reinterpret_cast<read_request *>(m_payload);
	}

	// Precondition: kind() == READ
	const read_request & get_read_request() const {
		return *reinterpret_cast<const read_request *>(m_payload);
	}

	// Precondition: kind() == WRITE
	write_request & get_write_request() {
		return *reinterpret_cast<write_request *>(m_payload);
	}

	// Precondition: kind() == WRITE
	const write_request & get_write_request() const {
		return *reinterpret_cast<const write_request *>(m_payload);
	}

	compressor_request_kind::type kind() const {
		return m_kind;
	}

private:
	void destruct() {
		switch (m_kind) {
			case compressor_request_kind::NONE:
				break;
			case compressor_request_kind::READ:
				get_read_request().~read_request();
				break;
			case compressor_request_kind::WRITE:
				get_write_request().~write_request();
				break;
		}
		m_kind = compressor_request_kind::NONE;
	}

	compressor_request_kind::type m_kind;

	char m_payload[sizeof(read_request) < sizeof(write_request) ? sizeof(write_request) : sizeof(read_request)];
};

class compressor_thread {
	class impl;
	impl * pimpl;

public:
	typedef boost::shared_ptr<compressor_buffer> buffer_t;
	typedef file_accessor::byte_stream_accessor<default_raw_file_accessor> file_accessor_t;
	typedef boost::mutex mutex_t;

	compressor_thread();
	~compressor_thread();

	mutex_t & mutex();

	// Locking: Caller must lock the thread (with a compressor_thread_lock).
	void request(compressor_request & r);

	void wait_for_request_done(compressor_thread_lock & l);

	void run();

	void stop(compressor_thread_lock & lock);
};

class compressor_thread_lock {
public:
	typedef boost::unique_lock<compressor_thread::mutex_t> lock_t;

	compressor_thread_lock(compressor_thread & c)
		: m_lock(c.mutex())
	{
	}

	lock_t & get_lock() {
		return m_lock;
	}

private:
	lock_t m_lock;
};

void read_request::wait(compressor_thread_lock & lock) {
	m_cond.wait(lock.get_lock());
}


}

#endif // TPIE_COMPRESSED_THREAD_H
