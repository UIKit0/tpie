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

#include <queue>
#include <snappy.h>
#include <tpie/compressed/thread.h>
#include <tpie/compressed/request.h>
#include <tpie/compressed/buffer.h>

namespace {

class block_header {
public:
	block_header()
		: m_payload(0)
	{
	}

	tpie::memory_size_type get_block_size() const {
		return static_cast<tpie::memory_size_type>(m_payload & BLOCK_SIZE_MASK);
	}

	// Precondition: 0 <= blockSize <= max_block_size()
	// Postcondition: get_block_size() == blockSize
	void set_block_size(tpie::memory_size_type blockSize) {
		m_payload &= ~BLOCK_SIZE_MASK;
		m_payload |= static_cast<tpie::uint32_t>(blockSize);
	}

	static tpie::memory_size_type max_block_size() {
		return BLOCK_SIZE_MAX;
	}

	bool operator==(const block_header & other) const {
		return m_payload == other.m_payload;
	}

	bool operator!=(const block_header & other) const { return !(*this == other); }

private:
	static const tpie::uint32_t BLOCK_SIZE_BITS = 24;
	static const tpie::uint32_t BLOCK_SIZE_MASK = (1 << BLOCK_SIZE_BITS) - 1;
	static const tpie::memory_size_type BLOCK_SIZE_MAX =
		static_cast<tpie::memory_size_type>(1 << BLOCK_SIZE_BITS) - 1;

	tpie::uint32_t m_payload;
};

}

namespace tpie {

/*static*/ stream_size_type compressor_thread::subtract_block_header(stream_size_type dataOffset) {
	return dataOffset - sizeof(block_header);
}

class compressor_thread::impl {
public:
	impl()
		: m_done(false)
	{
	}

	void stop(compressor_thread_lock & /*lock*/) {
		m_done = true;
		m_newRequest.notify_one();
	}

	bool request_valid(const compressor_request & r) {
		switch (r.kind()) {
			case compressor_request_kind::NONE:
				return false;
			case compressor_request_kind::READ:
			case compressor_request_kind::WRITE:
				return true;
		}
		throw exception("Unknown request type");
	}

	void run() {
		while (true) {
			compressor_thread_lock::lock_t lock(mutex());
			while (!m_done && m_requests.empty()) m_newRequest.wait(lock);
			if (m_done && m_requests.empty()) break;
			{
				compressor_request r = m_requests.front();
				m_requests.pop();
				lock.unlock();

				switch (r.kind()) {
					case compressor_request_kind::NONE:
						throw exception("Invalid request");
					case compressor_request_kind::READ:
						process_read_request(r.get_read_request());
						break;
					case compressor_request_kind::WRITE:
						process_write_request(r.get_write_request());
						break;
				}
			}
			lock.lock();
			m_requestDone.notify_all();
		}
	}

private:
	void process_read_request(read_request & rr) {
		const bool useCompression = rr.file_accessor().get_compressed();
		const bool backward = rr.read_direction() == direction::backward;
		if (backward && !useCompression) throw exception("backward && !useCompression");

		stream_size_type readOffset = rr.read_offset();
		if (!useCompression) {
			memory_size_type blockSize = rr.buffer()->size();
			if (blockSize > rr.buffer()->capacity()) {
				throw stream_exception("Internal error; blockSize > buffer capacity");
			}
			rr.file_accessor().read(readOffset, rr.buffer()->get(), blockSize);
			rr.buffer()->set_size(blockSize);
			compressor_thread_lock::lock_t lock(mutex());
			// Notify that reading has completed.
			rr.set_next_block_offset(1111111111111111111ull);
			return;
		}
		block_header blockHeader;
		block_header blockTrailer;
		memory_size_type blockSize;
		array<char> scratch;
		char * compressed;
		stream_size_type nextReadOffset;
		if (backward) {
			readOffset -= sizeof(blockTrailer);
			{
				memory_size_type nRead = rr.file_accessor().read(readOffset, &blockTrailer, sizeof(blockTrailer));
				if (nRead != sizeof(blockTrailer)) {
					throw exception("read failed to read right amount");
				}
			}
			blockSize = blockTrailer.get_block_size();
			if (blockSize == 0) {
				throw exception("Block size was unexpectedly zero");
			}
			scratch.resize(sizeof(blockHeader) + blockSize);
			readOffset -= scratch.size();
			{
				memory_size_type nRead = rr.file_accessor().read(readOffset, scratch.get(), scratch.size());
				if (nRead != scratch.size()) {
					throw exception("read failed to read right amount");
				}
			}
			compressed = scratch.get() + sizeof(blockHeader);
			memcpy(&blockHeader,
				   reinterpret_cast<block_header *>(scratch.get()),
				   sizeof(blockHeader));
			nextReadOffset = readOffset;
		} else {
			{
				memory_size_type nRead = rr.file_accessor().read(readOffset, &blockHeader, sizeof(blockHeader));
				if (nRead != sizeof(blockHeader)) {
					throw exception("read failed to read right amount");
				}
			}
			readOffset += sizeof(blockHeader);
			blockSize = blockHeader.get_block_size();
			if (blockSize == 0) {
				throw exception("Block size was unexpectedly zero");
			}
			scratch.resize(blockSize + sizeof(blockTrailer));
			{
				memory_size_type nRead = rr.file_accessor().read(readOffset, scratch.get(), scratch.size());
				if (nRead != scratch.size()) {
					throw exception("read failed to read right amount");
				}
			}
			compressed = scratch.get();
			memcpy(&blockTrailer,
				   reinterpret_cast<block_header *>(scratch.get() + scratch.size()) - 1,
				   sizeof(blockTrailer));
			nextReadOffset = readOffset + scratch.size();
		}
		if (blockHeader != blockTrailer) {
			throw exception("Block trailer is different from the block header");
		}

		size_t uncompressedLength;
		if (!snappy::GetUncompressedLength(compressed,
										   blockSize,
										   &uncompressedLength))
			throw stream_exception("Internal error; snappy::GetUncompressedLength failed");
		if (uncompressedLength > rr.buffer()->capacity()) {
			log_error() << uncompressedLength << ' ' << rr.buffer()->capacity() << std::endl;
			throw stream_exception("Internal error; snappy::GetUncompressedLength exceeds the block size");
		}
		rr.buffer()->set_size(uncompressedLength);
		snappy::RawUncompress(compressed,
							  blockSize,
							  reinterpret_cast<char *>(rr.buffer()->get()));

		compressor_thread_lock::lock_t lock(mutex());
		rr.set_next_block_offset(nextReadOffset);
	}

	void process_write_request(write_request & wr) {
		size_t inputLength = wr.buffer()->size();
		if (!wr.file_accessor().get_compressed()) {
			// Uncompressed case
			wr.file_accessor().write(wr.write_offset(), wr.buffer()->get(), wr.buffer()->size());
			return;
		}
		// Compressed case
		block_header blockHeader;
		block_header & blockTrailer = blockHeader;
		const memory_size_type maxBlockSize = snappy::MaxCompressedLength(inputLength);
		if (maxBlockSize > blockHeader.max_block_size())
			throw exception("process_write_request: MaxCompressedLength > max_block_size");
		array<char> scratch(sizeof(blockHeader) + maxBlockSize + sizeof(blockTrailer));
		memory_size_type blockSize;
		snappy::RawCompress(reinterpret_cast<const char *>(wr.buffer()->get()),
							inputLength,
							scratch.get() + sizeof(blockHeader),
							&blockSize);
		blockHeader.set_block_size(blockSize);
		memcpy(scratch.get(), &blockHeader, sizeof(blockHeader));
		memcpy(scratch.get() + sizeof(blockHeader) + blockSize, &blockTrailer, sizeof(blockTrailer));
		const memory_size_type writeSize = sizeof(blockHeader) + blockSize + sizeof(blockTrailer);
		if (!wr.should_append()) {
			log_debug() << "Truncate to " << wr.write_offset() << std::endl;
			wr.file_accessor().truncate_bytes(wr.write_offset());
			log_debug() << "File size is now " << wr.file_accessor().file_size() << std::endl;
		}
		{
			compressor_thread_lock::lock_t lock(mutex());
			wr.set_block_info(wr.file_accessor().file_size(), writeSize);
		}
		wr.file_accessor().append(scratch.get(), writeSize);
	}

public:
	mutex_t & mutex() {
		return m_mutex;
	}

	void request(const compressor_request & r) {
		if (!request_valid(r))
			throw exception("Invalid request");

		m_requests.push(r);
		m_requests.back().get_request_base().initiate_request();
		m_newRequest.notify_one();
	}

	void wait_for_request_done(compressor_thread_lock & l) {
		m_requestDone.wait(l.get_lock());
	}

private:
	mutex_t m_mutex;
	std::queue<compressor_request> m_requests;
	boost::condition_variable m_newRequest;
	boost::condition_variable m_requestDone;
	bool m_done;
};

} // namespace tpie

namespace {

tpie::compressor_thread the_compressor_thread;
boost::thread the_compressor_thread_handle;
bool compressor_thread_already_finished = false;

void run_the_compressor_thread() {
	the_compressor_thread.run();
}

} // unnamed namespace

namespace tpie {

compressor_thread & the_compressor_thread() {
	return ::the_compressor_thread;
}

void init_compressor() {
	if (the_compressor_thread_handle.get_id() != boost::thread::id()) {
		log_debug() << "Attempted to initiate compressor thread twice" << std::endl;
		return;
	}
	boost::thread t(run_the_compressor_thread);
	the_compressor_thread_handle.swap(t);
	compressor_thread_already_finished = false;
}

void finish_compressor() {
	if (the_compressor_thread_handle.get_id() == boost::thread::id()) {
		if (compressor_thread_already_finished) {
			log_debug() << "Compressor thread already finished" << std::endl;
		} else {
			log_debug() << "Attempted to finish compressor thread that was never initiated" << std::endl;
		}
		return;
	}
	{
		compressor_thread_lock lock(the_compressor_thread());
		the_compressor_thread().stop(lock);
	}
	the_compressor_thread_handle.join();
	boost::thread t;
	the_compressor_thread_handle.swap(t);
	compressor_thread_already_finished = true;
}

compressor_thread::compressor_thread()
	: pimpl(new impl)
{
}

compressor_thread::~compressor_thread() {
	delete pimpl;
}

compressor_thread::mutex_t & compressor_thread::mutex() {
	return pimpl->mutex();
}

void compressor_thread::request(compressor_request & r) {
	pimpl->request(r);
}

void compressor_thread::run() {
	pimpl->run();
}

void compressor_thread::wait_for_request_done(compressor_thread_lock & l) {
	pimpl->wait_for_request_done(l);
}

void compressor_thread::stop(compressor_thread_lock & lock) {
	pimpl->stop(lock);
}

}
