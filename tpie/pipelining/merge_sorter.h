// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :
// Copyright 2012, 2014, The TPIE development team
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
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TPIE. If not, see <http://www.gnu.org/licenses/>

#ifndef __TPIE_PIPELINING_MERGE_SORTER_H__
#define __TPIE_PIPELINING_MERGE_SORTER_H__

#include <tpie/tpie.h>
#include <tpie/pipelining/sort_parameters.h>
#include <tpie/file_stream.h>
#include <tpie/tempname.h>
#include <tpie/blocking_queue.h>
#include <tpie/parallel_sort.h>
#include <tpie/internal_priority_queue.h>
#include <tpie/internal_vector.h>
#include <tpie/tournament_tree.h>
#include <deque>
#include <boost/atomic.hpp>

namespace tpie {

///////////////////////////////////////////////////////////////////////////////
/// Merge sorting consists of three phases.
///
/// 1. Sorting and forming runs
/// 2. Merging runs
/// 3. Final merge and report
///
/// If the number of elements received during phase 1 is less than the length
/// of a single run, we are in "report internal" mode, meaning we do not write
/// anything to disk. This causes phase 2 to be a no-op and phase 3 to be a
/// simple array traversal.
///////////////////////////////////////////////////////////////////////////////
template <typename T, bool UseProgress, typename pred_t = std::less<T> >
class merge_sorter {
private:
	typedef internal_vector<T> run_container_type;

	struct block { // used in the producer/consumer-pattern in phase 2 and 3
		struct mode { // C++03-style enum
			enum type {
				data,
				run_signal,
				terminate_signal
			};
		};

		block(typename mode::type m) : m_mode(m) {
			if(m_mode == block::mode::data) {
				m_data = tpie_new<run_container_type>(get_block_size() / sizeof(T));
			}
		}

		~block() {
			if(m_mode == block::mode::data) {
				tpie_delete(m_data);
			}
		}

		static memory_size_type memory_usage() {
			return sizeof(block)
			+ (get_block_size() / sizeof(T)) * sizeof(T);
		}

		typename mode::type m_mode;
		bool last_block;
		memory_size_type file;
		run_container_type * m_data; // If mode is set to data, this pointer will be set
	};

	struct tournament_leaf {
		block * block_pointer;
		bool last_block;
		T * begin;
		T * end;

		T smallest_element;
	};

	typedef std::vector<tournament_leaf> leaves;

	template<typename S>
	class predwrap {
	public:
		typedef const std::pair<T, S>& item_type;
		typedef item_type first_argument_type;
		typedef item_type second_argument_type;
		typedef bool result_type;

		predwrap(pred_t & pred) : m_pred(pred) {}

		bool operator()(item_type a, item_type b) {
			return m_pred(a.first, b.first);
		}
	private:
		pred_t m_pred;
	};

	class tourn_pred {
	public:
		typedef const tournament_leaf& item_type;
		typedef item_type first_argument_type;
		typedef item_type second_argument_type;
		typedef bool result_type;

		tourn_pred(pred_t & pred) : m_pred(pred) {}

		tourn_pred(const tourn_pred & other) : m_pred(other.m_pred) {}

		bool operator()(const tournament_leaf& a, const tournament_leaf& b) {
			return m_pred(a.smallest_element, b.smallest_element);
		}
	private:
		pred_t m_pred;
	};

	enum reporting_mode {
		REPORTING_MODE_INTERNAL, // do not write anything to disk. Phase 2 will be a no-op phase and phase 3 is simple array traversal
		REPORTING_MODE_EXTERNAL
	};

	enum state_type {
		STATE_PARAMETERS,
		STATE_RUN_FORMATION,
		STATE_MERGE,
		STATE_REPORT
	};
public:
	typedef boost::shared_ptr<merge_sorter> ptr;
	typedef progress_types<UseProgress> Progress;

	static const memory_size_type maximumFanout = 250; // This is the max number of runs to merge at a time when running a k-way merge.
	static const memory_size_type bufferCount = 2; // This is the number of buffers to be used during phase 1
	static const memory_size_type writebufferCount = 3; // the number of write buffers used in phase 2

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Create a new merger_sorter object with the given predicate
	///////////////////////////////////////////////////////////////////////////////
	merge_sorter(pred_t pred = pred_t())
	: m_parametersSet(false)
	, m_pred(pred)
	, m_reporting_mode(REPORTING_MODE_EXTERNAL)
	, m_state(STATE_PARAMETERS)
	, m_runsPushed(0)
	, m_itemsLeft(0)
	, m_itemsPushed(0)
	, m_largestElementSet(false)
		, m_tournamentTree(tourn_pred(m_pred)) // set the size to 0 for now
	{
		m_parameters.memoryPhase1 = 0;
		m_parameters.memoryPhase2 = 0;
		m_parameters.memoryPhase3 = 0;
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Enable setting run length and fanout manually (for testing
	/// purposes).
	///////////////////////////////////////////////////////////////////////////
	void set_parameters(memory_size_type runLength, memory_size_type fanout) {
		tp_assert(m_state == STATE_PARAMETERS, "Merge sorting already begun");
		m_parameters.runLength = runLength;
		m_parameters.internalReportThreshold = runLength/3;
		m_parameters.fanout = m_parameters.finalFanout = fanout;

		m_parametersSet = true;
	}

private:
	static memory_size_type phase2_fanout_memory_usage(memory_size_type fanout) {
		return sizeof(merge_sorter<T, UseProgress, pred_t>)
		+ fanout * block::memory_usage();
	}

	static memory_size_type calculate_phase2_fanout(memory_size_type memory) {
		// solves the equation given in phase2_fanout_memory_usage
		return (memory - sizeof(merge_sorter<T, UseProgress, pred_t>)) / block::memory_usage();
	}


	static memory_size_type phase3_fanout_memory_usage(memory_size_type fanout) {
		memory_size_type overhead = sizeof(merge_sorter<T, UseProgress, pred_t>) + tpie::internal_priority_queue<std::pair<T, block*>, predwrap<block*> >::memory_overhead();
		memory_size_type coefficient = block::memory_usage() + tpie::internal_priority_queue<std::pair<T, block*>, predwrap<block*> >::memory_coefficient();

		return overhead + fanout * coefficient;
	}

	static memory_size_type calculate_phase3_fanout(memory_size_type memory) {
		// solves the equation given in phase3_fanout_memory_usage
		memory_size_type overhead = sizeof(merge_sorter<T, UseProgress, pred_t>) + tpie::internal_priority_queue<std::pair<T, block*>, predwrap<block*> >::memory_overhead();
		memory_size_type coefficient = block::memory_usage() + tpie::internal_priority_queue<std::pair<T, block*>, predwrap<block*> >::memory_coefficient();

		return (memory - overhead) / coefficient;
	}

	void calculate_parameters() {
		tp_assert(m_state == STATE_PARAMETERS, "Merge sorting already begun");

		///////////////////////////////////////////////////////////////////////////////
		/// Phase 2
		/// The fanout is determined by the size of the merge heap and the stream
		/// memory usage
		///////////////////////////////////////////////////////////////////////////////

		memory_size_type availableReadMemory = m_parameters.memoryPhase2; // the memory available for read buffers in phase 2
		if(m_parameters.memoryPhase2 >= block::memory_usage() * writebufferCount)
			availableReadMemory -= block::memory_usage() * writebufferCount;
		else
			log_debug() << "Not enough memory for writeBuffers in phase 1" << block::memory_usage() * writebufferCount << " > " << m_parameters.memoryPhase2;

		m_parameters.fanout = calculate_phase2_fanout(availableReadMemory); // the write buffers used in phase 2 should be taken into account. These are not used in phase 3.
		if(phase2_fanout_memory_usage(m_parameters.fanout) > availableReadMemory) {
			log_debug() << "Not enough memory for fanout " << m_parameters.fanout << " "
						<< phase2_fanout_memory_usage(m_parameters.fanout) << " > "
						<< availableReadMemory << std::endl;

			m_parameters.memoryPhase2 = phase2_fanout_memory_usage(m_parameters.fanout);
		}


		///////////////////////////////////////////////////////////////////////////////
		/// Phase 1
		/// The run length is determined by the number of items we can hold in memory
		/// and the number of buffers we are using
		///////////////////////////////////////////////////////////////////////////////

		m_parameters.runLength = (m_parameters.memoryPhase1 - sizeof(merge_sorter<T, UseProgress, pred_t>)) / (sizeof(T) * bufferCount);
		m_desiredSize = m_parameters.runLength;

		// if we receive less items than internalReportThreshold, internal report mode will be used(no I/O)
		m_parameters.internalReportThreshold = std::min(m_parameters.memoryPhase1, std::min(m_parameters.memoryPhase2, m_parameters.memoryPhase3)) / sizeof(T);
		if(m_parameters.internalReportThreshold > m_parameters.runLength/3)
			m_parameters.internalReportThreshold = m_parameters.runLength/3;

		///////////////////////////////////////////////////////////////////////////////
		/// Phase 3
		///////////////////////////////////////////////////////////////////////////////
		m_parameters.finalFanout = calculate_phase3_fanout(m_parameters.memoryPhase3); // there are no write buffers used in phase 3, only read buffers.

		if(phase3_fanout_memory_usage(m_parameters.finalFanout) > m_parameters.memoryPhase3) {
			log_debug() << "Not enough memory for fanout "
						<< m_parameters.finalFanout
						<< "! (" << m_parameters.memoryPhase3
						<< " < " << phase3_fanout_memory_usage(m_parameters.finalFanout) << ")\n";

			m_parameters.memoryPhase3 = phase3_fanout_memory_usage(m_parameters.finalFanout);
		}

		///////////////////////////////////////////////////////////////////////////////
		/// Final
		///////////////////////////////////////////////////////////////////////////////
		m_parametersSet = true;
	}
public:

	///////////////////////////////////////////////////////////////////////////
	/// \brief Calculate parameters from given memory amount.
	/// \param m Memory available for phase 2, 3 and 4
	///////////////////////////////////////////////////////////////////////////
	void set_available_memory(memory_size_type m) {
		m_parameters.memoryPhase1 = m_parameters.memoryPhase2 = m_parameters.memoryPhase3 = m;
		calculate_parameters();
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Calculate parameters from given memory amount.
	/// \param m1 Memory available for phase 1
	/// \param m2 Memory available for phase 2
	/// \param m3 Memory available for phase 3
	///////////////////////////////////////////////////////////////////////////
	void set_available_memory(memory_size_type m1, memory_size_type m2, memory_size_type m3) {
		m_parameters.memoryPhase1 = m1;
		m_parameters.memoryPhase2 = m2;
		m_parameters.memoryPhase3 = m3;

		calculate_parameters();
	}

private:
	///////////////////////////////////////////////////////////////////////////////
	/// \brief Helper for set_phase_?_memory
	///////////////////////////////////////////////////////////////////////////////
	void maybe_calculate_parameters() {
		if(m_parameters.memoryPhase1 > 0 && m_parameters.memoryPhase2 > 0 && m_parameters.memoryPhase3 > 0) {
			calculate_parameters();
		}
	}
public:
	///////////////////////////////////////////////////////////////////////////////
	/// \brief Calculate parameters from given memory amount for phase 1.
	/// \param mem Memory avaiable for phase 1
	///////////////////////////////////////////////////////////////////////////////
	void set_phase_1_memory(memory_size_type mem) {
		m_parameters.memoryPhase1 = mem;
		maybe_calculate_parameters();
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Calculate parameters from given memory amount for phase 2.
	/// \param mem Memory avaiable for phase 2
	///////////////////////////////////////////////////////////////////////////////
	void set_phase_2_memory(memory_size_type mem) {
		m_parameters.memoryPhase2 = mem;
		maybe_calculate_parameters();
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Calculate parameters from given memory amount for phase 3.
	/// \param mem Memory avaiable for phase 3
	///////////////////////////////////////////////////////////////////////////////
	void set_phase_3_memory(memory_size_type mem) {
		m_parameters.memoryPhase3 = mem;
		maybe_calculate_parameters();
	}

	// This thread is only used during phase 1 for sorting runs
	void phase1_sort_thread() {
		// initialization
		for(memory_size_type i = 1; i < bufferCount; ++i)
			m_emptyBuffers.push(tpie_new<run_container_type>(m_parameters.runLength));

		m_IOThread = boost::thread(boost::bind(&merge_sorter::io_thread, this));

		// Sorting
		while(true) {
			run_container_type * run = m_fullBuffers.pop();
			if(run == NULL) {
				m_sortedBuffers.push(NULL);
				break;
			}

			tpie::parallel_sort(run->begin(), run->end(), m_pred);

			m_sortedBuffers.push(run);
		}
	}

	// this thread is used during all 3 fases to perform IO operations
	void io_thread() {
		// phase 1
		while(true) {
			run_container_type * run = m_sortedBuffers.pop();
			if(run == NULL) {
				break;
			}

			if(!m_largestElementSet || m_pred(m_largestElement, run->back())) {
				m_largestElement = run->back();
				m_largestElementSet = true;
			}

			temp_file runFile;
			file_accessor::raw_file_accessor out;
			out.open_wo(runFile.path()); // open for writing
			out.write_i((void*) (&run->front() + 1), (run->size() - 1) * sizeof(T));
			out.close_i();

			m_runFiles.push_back(runFile);
			m_smallestElements.push_back(run->front());
			m_runFileSizes.push_back(run->size() - 1);
			run->clear();
			m_emptyBuffers.push(run); // push a new empty buffer
		}

		// phase 2

		// phase 3
		tpie::array<file_accessor::raw_file_accessor> in(m_runFiles.size());
		memory_size_type m = 0; // the number of elements left to read

		for(memory_size_type i = 0; i < m_runFiles.size(); ++i) {
			in[i].open_ro(m_runFiles[i].path());
			m += m_runFileSizes[i];
		}

		memory_size_type n = 2;
		while(n < m_runFiles.size()) {
			n *= 2;
		}

		tpie::array<T> runs(n);
		std::copy(m_smallestElements.begin(), m_smallestElements.end(), runs.begin());

		for(memory_size_type i = m_runFiles.size(); i < n; ++i) {
			runs[i] = m_largestElement;
		}

		bits::tournament_tree<typename tpie::array<T>::iterator, pred_t> tree(runs.begin(), runs.end(), m_pred);

		while(m > 0) {
			block * b;
			{
				boost::mutex::scoped_lock lock(m_mutex);
				while(m_emptyReadBuffers.empty()) {
					m_condition.wait(lock);
				}

				b = m_emptyReadBuffers.front();
				m_emptyReadBuffers.pop();
			}

			memory_size_type run = tree.top();
			memory_size_type size = std::min(get_block_size() / sizeof(T), m_runFileSizes[run]);

			if(size < get_block_size() / sizeof(T)) {
				// this is the last block of the run. The last element therefore need to be included
				b->m_data->resize(size+1);
				b->m_data->front() = runs[run];
				in[run].read_i((void*) ((&b->m_data->front())+1), size * sizeof(T));

				runs[run] = m_largestElement;

				m_runFileSizes[run] -= size;
				m -= size;

				b->last_block = true;
			}
			else {
				b->m_data->resize(size); // the extra element is only used temporarily in this thread
				b->m_data->front() = runs[run]; // the smallest element
				in[run].read_i((void*) ((&b->m_data->front())+1), (size-1) * sizeof(T));

				runs[run] = b->m_data->back();
				b->m_data->pop_back();

				m_runFileSizes[run] -= size-1;
				m -= size-1;

				b->last_block = false;
			}

			tree.update_key(run);
			b->file = run;

			boost::mutex::scoped_lock lock(m_mutex);
			m_fullReadBuffers.push(b);
			lock.unlock(); // avoid problems if a waiting thread immediately wakes up
			m_condition.notify_one();
		}
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Initiate phase 1: Formation of input runs.
	///////////////////////////////////////////////////////////////////////////
	void begin() {
		log_debug() << "Starting merge sorter phase 1" << std::endl;
		tp_assert(m_parametersSet, "Parameters have not been set");
		m_state = STATE_RUN_FORMATION;

		m_currentRun = tpie_new<run_container_type>(m_parameters.runLength);
		m_sortThread = boost::thread(boost::bind(&merge_sorter::phase1_sort_thread, this)); // perform the rest of the initialization in the sort thread
	}

private:
	void push_run() {
		m_fullBuffers.push(m_currentRun);
		m_currentRun = m_emptyBuffers.pop();
		++m_runsPushed;
		if(m_runsPushed == 1)
			m_desiredSize = m_parameters.runLength/2;
		else
			m_desiredSize = m_parameters.runLength;
	}

public:
	///////////////////////////////////////////////////////////////////////////
	/// \brief Push item to merge sorter during phase 1.
	///////////////////////////////////////////////////////////////////////////
	void push(const T & item) {
		if(m_currentRun->size() == m_desiredSize) push_run(); // push the run to the sort thread and fetch a new empty buffer

		m_currentRun->push_back(item);
		++m_itemsPushed;
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief End phase 1.
	///////////////////////////////////////////////////////////////////////////
	void end() {
		// Set to internal reporting mode if possible

		m_itemsLeft = m_itemsPushed;

		if(m_runsPushed == 0 && m_itemsPushed <= m_parameters.internalReportThreshold) {
			/*log_info() << "Internal reporting mode set. " << std::endl;
			log_info() << "Capacity: " << m_emptyBuffers.front()->capacity() << std::endl;
			log_info() << "Size: " << m_emptyBuffers.front()->size() << std::endl;
			log_info() << "Threshold: " << m_parameters.internalReportThreshold << std::endl;
			log_info() << "Items pushed: " << m_itemsPushed << std::endl;*/

			// no buffer has been queued for sorting or writing yet and the number size is small enough.
			// use internal report mode

			m_fullBuffers.push(NULL);

			if(!m_currentRun->empty())
				tpie::parallel_sort(m_currentRun->begin(), m_currentRun->end(), m_pred);
			else {
				m_emptyBuffers.push(m_currentRun);
				m_currentRun = NULL;
			}

			m_reporting_mode = REPORTING_MODE_INTERNAL;
		}
		else {
			/*log_info() << "External reporting mode set. " << std::endl;
			log_info() << "Capacity: " << m_emptyBuffers.front()->capacity() << std::endl;
			log_info() << "Size: " << m_emptyBuffers.front()->size() << std::endl;
			log_info() << "Threshold: " << m_parameters.internalReportThreshold << std::endl;
			log_info() << "Items pushed: " << m_itemsPushed << std::endl;*/

			// use external report mode

			if(!m_currentRun->empty())
				m_fullBuffers.push(m_currentRun);
			else
				m_emptyBuffers.push(m_currentRun);

			m_currentRun = NULL;
			m_fullBuffers.push(NULL);
		}


		m_state = STATE_MERGE;

		m_sortThread.join();

		for(memory_size_type i = (m_currentRun != NULL ? 1 : 0); i < bufferCount; ++i) {
			tpie_delete(m_emptyBuffers.pop());
		}

		tp_assert(m_emptyBuffers.empty(), "A wild write buffer has appeared.");

		if(m_itemsPushed == 0) {
			boost::mutex::scoped_lock lock(m_mutex);
			m_fullWriteBuffers.push(tpie_new<block>(block::mode::terminate_signal)); // push a run to signal thread termination
			lock.unlock(); // avoid problems if a waiting thread immediately wakes up
			m_condition.notify_one();
			m_IOThread.join();
		}
	}

private:
	///////////////////////////////////////////////////////////////////////////////
	/// \brief Merges the first n runs in the run deque and pushes a new runfile
	/// to the end
	///////////////////////////////////////////////////////////////////////////////
	void merge_runs(memory_size_type n, typename Progress::base & pi) {
		/*log_debug() << "Merging " << n << std::endl;
		if(n > m_runFiles.size()) {
			n = m_runFiles.size();
		}

		std::vector<std::pair<T, memory_size_type> > sigma;
		for(memory_size_type i = 0; i < n; ++i) {
			for(typename std::vector<T>::iterator j = m_phi[i].begin(); j < m_phi[i].end(); ++j) {
				sigma.push_back(std::make_pair(*j, i));
			}
		}

		tpie::parallel_sort(sigma.begin(), sigma.end(), predwrap<memory_size_type>(m_pred));

		for(typename std::vector<std::pair<T, memory_size_type> >::iterator i = sigma.begin(); i != sigma.end(); ++i) {
			m_readJobs.push(i->second);
		}

		memory_size_type next_block = 0;
		tpie::internal_priority_queue<std::pair<T, block*>, predwrap<block*> > queue(m_parameters.fanout, predwrap<block*>(m_pred));

		block * write_buffer = m_emptyWriteBuffers.pop();

		while(true) {
			if(next_block == sigma.size()) // all blocks have been processed
				break;

			do {
				block * b = m_fullReadBuffers.pop();
				pi.step();
				queue.push(std::make_pair(b->m_data->front(), b));
				b->m_data->pop_front();
				++next_block;
			} while(!m_fullReadBuffers.empty());

			while(!queue.empty() && !(next_block < sigma.size() && m_pred(sigma[next_block].first, queue.top().first))) {
				const std::pair<T, block*> & top = queue.top();
				write_buffer->m_data->push_back(top.first);
				if(write_buffer->m_data->size() == get_block_size()/sizeof(T)) { // replace the full write buffer
					m_fullWriteBuffers.push(write_buffer);
					write_buffer = m_emptyWriteBuffers.pop();
				}

				block * b = top.second;
				if(b->m_data->empty()) {
					b->m_data->clear();
					m_emptyReadBuffers.push(b);
					queue.pop();
				}
				else {
					queue.pop_and_push(std::make_pair(b->m_data->front(), b));
					b->m_data->pop_front();
				}
			}
		}

		if(!write_buffer->m_data->empty()) {
			m_fullWriteBuffers.push(write_buffer);
		}
		else {
			m_emptyWriteBuffers.push(write_buffer);
		}

		m_fullWriteBuffers.push(tpie_new<block>(block::mode::run_signal));

		while(m_queuedWriteJobs > 0)
			boost::this_thread::yield();

		for(memory_size_type i = 0; i < n; ++i) {
			m_runFiles.pop_front();
			m_runFileSizes.pop_front();
			m_phi.pop_front();
		}*/
	}

public:

	///////////////////////////////////////////////////////////////////////////
	/// \brief Perform phase 2: Performing all merges in the merge tree
	///////////////////////////////////////////////////////////////////////////
	void calc(typename Progress::base & pi) {
		log_debug() << "Starting merge sorter phase 2" << std::endl;
		if(m_reporting_mode == REPORTING_MODE_INTERNAL) {
			pi.init(1);
			pi.step();
			pi.done();
		}
		else {
			if(m_runFiles.size() > m_parameters.finalFanout) {
				tp_assert(false, "Phase 2 is not implemented yet");
			/*	memory_size_type treeHeight = static_cast<int> (
					ceil(log(static_cast<float>(m_runFiles.size()))) / ceil(log(static_cast<float>(m_parameters.fanout)))
				); // the number of merge 'rounds' to be performed.

				pi.init(treeHeight * m_itemsPushed / get_block_size());

				// allocate buffers
				for(memory_size_type i = 0; i < m_parameters.fanout; ++i)
					m_emptyReadBuffers.push(tpie_new<block>(block::mode::data));

				for(memory_size_type i = 0; i < writebufferCount; ++i)
					m_emptyWriteBuffers.push(tpie_new<block>(block::mode::data));

				// perform merges
				do {
					merge_runs(m_parameters.fanout, pi);
				} while(m_runFiles.size() > m_parameters.finalFanout);

				// destroy the buffers
				for(memory_size_type i = 0; i < m_parameters.fanout; ++i)
					tpie_delete(m_emptyReadBuffers.pop());

				for(memory_size_type i = 0; i < writebufferCount; ++i)
					tpie_delete(m_emptyWriteBuffers.pop());

				pi.done();*/
			}
			else {
				pi.init(1);
				pi.step();
				pi.done();
			}

		}

		m_state = STATE_REPORT;
		//tpie::log_info() << "Finished phase 2" << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief If in internal reporting mode: write data to disk. The allocated
	/// memory is needed elsewhere
	///////////////////////////////////////////////////////////////////////////////
	void evacuate() {
		tp_assert(m_state == STATE_MERGE || m_state == STATE_REPORT, "Wrong phase");

		if(m_reporting_mode == REPORTING_MODE_INTERNAL) {
			if(m_currentRun != NULL) { // the buffer is non-empty -> write the buffer to disk
				/*m_reporting_mode = REPORTING_MODE_EXTERNAL; // write the buffer to disk and use external reporting mode
				temp_file runFile;
				file_accessor::raw_file_accessor out;
				out.open_wo(runFile.path()); // open for writing
				out.write_i((void*) &m_currentRun->front(), m_currentRun->size() * sizeof(T));
				out.close_i();

				m_runFiles.push_back(runFile);
				m_runFileSizes.push_back(m_currentRun->size() * sizeof(T));

				m_phi.push_back(std::vector<T>());
				std::vector<T> & phi = m_phi.back();
				for(memory_size_type i = 0; i < m_currentRun->size(); i += get_block_size()/sizeof(T)) phi.push_back((*m_currentRun)[i]);

				tpie_delete(m_currentRun);*/
				tp_assert(false, "not implemented yet");
			}
		}
		else {
			// nothing to do in external reporting mode
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Evacuate if we are in the merging state
	///////////////////////////////////////////////////////////////////////////////
	void evacuate_before_merging() {
		if(m_state == STATE_MERGE)
		evacuate();
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Evacuate if we are in the reporting state
	///////////////////////////////////////////////////////////////////////////////
	void evacuate_before_reporting() {
		if(m_state == STATE_REPORT && (m_reporting_mode == REPORTING_MODE_INTERNAL || m_itemsLeft == m_itemsPushed))
			evacuate();
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief In phase 3, return true if there are more items in the final merge
	/// phase.
	///////////////////////////////////////////////////////////////////////////
	bool can_pull() const {
		tp_assert(m_state == STATE_REPORT, "Wrong phase");
		return m_itemsLeft != 0;
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Begin phase 3
	///////////////////////////////////////////////////////////////////////////////
	void pull_begin() {
		log_debug() << "Starting merge sorter phase 3" << std::endl;
		if(m_reporting_mode == REPORTING_MODE_INTERNAL) { // nothing to do in internal reporting mode
			return;
		}

		// initialize tournament tree
		// Let the size be the next power of 2
		memory_size_type n = 2;
		while(n < m_runFiles.size()) {
			n *= 2;
		}
		m_leaves.resize(n);

		for(memory_size_type i = 0; i < m_runFiles.size(); ++i) {
			tournament_leaf & leaf = m_leaves[i];
			leaf.smallest_element = m_smallestElements[i];
			leaf.last_block	= false;
			leaf.begin = leaf.end = NULL;
			leaf.block_pointer = NULL;
		}

		for(memory_size_type i = m_runFiles.size(); i < n; ++i) {
			tournament_leaf & leaf = m_leaves[i];
			leaf.smallest_element = m_largestElement;
			leaf.last_block	= true;
			leaf.begin = leaf.end = NULL;
			leaf.block_pointer = NULL;
		}

		m_tournamentTree = bits::tournament_tree<typename leaves::iterator, tourn_pred>(m_leaves.begin(), m_leaves.end(), tourn_pred(m_pred));

		// initialize IO
		{
			boost::mutex::scoped_lock lock(m_mutex);
			for(memory_size_type i = 0; i < m_parameters.finalFanout; ++i)
				m_emptyReadBuffers.push(tpie_new<block>(block::mode::data));
			lock.unlock(); // avoid problems if a waiting thread immediately wakes up
			m_condition.notify_one();
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief End phase 3
	///////////////////////////////////////////////////////////////////////////////
	void pull_end() {
		{
			boost::mutex::scoped_lock lock(m_mutex);
			m_fullWriteBuffers.push(tpie_new<block>(block::mode::terminate_signal)); // push a run to signal thread termination
			lock.unlock(); // avoid problems if a waiting thread immediately wakes up
			m_condition.notify_one();
		}

		m_IOThread.join();

		if(m_reporting_mode == REPORTING_MODE_INTERNAL)  {
			tpie_delete(m_currentRun);
			return;
		}

		// it is no longer neccesary to acquire the lock, since the IO thread has been stopped
		for(memory_size_type i = 0; i < m_parameters.finalFanout; ++i) {
			if(m_leaves[i].block_pointer != NULL) {
				m_emptyReadBuffers.push(m_leaves[i].block_pointer);
			}
		}

		for(memory_size_type i = 0; i < m_parameters.finalFanout; ++i) {
			tpie_delete(m_emptyReadBuffers.front());
			m_emptyReadBuffers.pop();
		}

		m_runFiles.clear();
	}

	void fetch_blocks(memory_size_type i) {
		tournament_leaf	& leaf = m_leaves[i];
		if(leaf.last_block) {
			if(leaf.block_pointer != NULL) {
				{
					boost::mutex::scoped_lock lock(m_mutex);
					m_emptyReadBuffers.push(leaf.block_pointer);
					lock.unlock(); // avoid problems if a waiting thread immediately wakes up
					m_condition.notify_one();
				}

				leaf.block_pointer = NULL;
				leaf.end = leaf.begin = &m_largestElement;
				++leaf.end;
			}

			return;
		}

		boost::mutex::scoped_lock lock(m_mutex);
		while(true) {
			while(m_fullReadBuffers.empty()) {
				m_condition.wait(lock);
			}
			block * b = m_fullReadBuffers.front(); m_fullReadBuffers.pop();
			tournament_leaf & leaf = m_leaves[b->file];
			tp_assert(leaf.begin == leaf.end, "Not empty");
			tp_assert(b != NULL, "the block is null");
			tp_assert(b->m_data != NULL, "the data is null");
			tp_assert(!b->m_data->empty(), "block is empty");
			leaf.begin = &((*b->m_data)[0]);
			leaf.last_block = b->last_block;
			leaf.end = &((*b->m_data)[0]) + b->m_data->size();
			if(leaf.block_pointer != NULL) {
				m_emptyReadBuffers.push(leaf.block_pointer);
				m_condition.notify_one();
			}
			leaf.block_pointer = b;

			if(b->file == i)
				break;
		}
	}

	T pull() {
		if(m_reporting_mode == REPORTING_MODE_INTERNAL) {
			return (*m_currentRun)[m_itemsPushed-(m_itemsLeft--)];
		}

		memory_size_type index = m_tournamentTree.top();
		tournament_leaf & leaf = m_leaves[index];

		T res = leaf.smallest_element;

		if(leaf.begin == leaf.end) {
			fetch_blocks(index);
		}

		leaf.smallest_element = *leaf.begin;
		++leaf.begin;
		m_tournamentTree.update_key(index);
		--m_itemsLeft;

		return res;
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief Return the number of items in the sorter
	///////////////////////////////////////////////////////////////////////////////
	stream_size_type item_count() {
		return m_itemsLeft;
	}
public:
	static memory_size_type memory_usage_phase_1(const sort_parameters & params) {
		return sizeof(merge_sorter<T, UseProgress, pred_t>) + params.runLength * sizeof(T) * bufferCount;
	}

	static memory_size_type minimum_memory_phase_1() {
		// a runlength of 1
		return sizeof(T) * bufferCount + sizeof(merge_sorter<T, UseProgress, pred_t>);
	}

	static memory_size_type memory_usage_phase_2(const sort_parameters & params) {
		return phase2_fanout_memory_usage(params.fanout) + writebufferCount * get_block_size();
	}

	static memory_size_type minimum_memory_phase_2() {
		return phase2_fanout_memory_usage(0) + writebufferCount * get_block_size();
	}

	static memory_size_type memory_usage_phase_3(const sort_parameters & params) {
		return phase3_fanout_memory_usage(params.finalFanout);
	}

	static memory_size_type minimum_memory_phase_3() {
		// The minimum amount of memory used is when the fanout is zero
		return phase3_fanout_memory_usage(0);
	}

	static memory_size_type maximum_memory_phase_3() {
		// The maximum amount of memory is used when the fanout is the maximum possible fanout
		return phase3_fanout_memory_usage(maximumFanout);
	}

	///////////////////////////////////////////////////////////////////////////////
	/// \brief The memory usage when the sorter is evacuated.
	///////////////////////////////////////////////////////////////////////////////
	memory_size_type evacuated_memory_usage() const {
		return sizeof(merge_sorter<T, UseProgress, pred_t>);
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Set upper bound on number of items pushed.
	///
	/// If the number of items to push is less than the size of a single run,
	/// this method will decrease the run size to that.
	/// This may make it easier for the sorter to go into internal reporting
	/// mode.
	///////////////////////////////////////////////////////////////////////////
	void set_items(stream_size_type) {
		// TODO: Avoid going into external report mode if possible
	}

private:
	sort_parameters m_parameters;
	bool m_parametersSet;
	pred_t m_pred;
	reporting_mode m_reporting_mode;
	state_type m_state;
	memory_size_type m_runsPushed;
	memory_size_type m_itemsLeft;
	memory_size_type m_itemsPushed;
	memory_size_type m_desiredSize;

	run_container_type * m_currentRun;
	T m_largestElement;
	bool m_largestElementSet;
	boost::thread m_IOThread; // The thread in phase 1 used to write run formations to file

	// as far as i can see these do not need to be locked with a mutex
	std::deque<temp_file> m_runFiles;
	std::deque<memory_size_type> m_runFileSizes;

	// phase 1 specific
	bits::blocking_queue<run_container_type *> m_emptyBuffers; // the buffers that are to be consumed by the push method
	bits::blocking_queue<run_container_type *> m_fullBuffers; // the buffers that are to be consumed by the sorting thread
	bits::blocking_queue<run_container_type *> m_sortedBuffers; // the buffers that are to be consumed by the write thread
	boost::thread m_sortThread; // The thread in phase 1 used to sort run formations

	// phase 2 + 3 specific
	std::deque<T> m_smallestElements;
	std::queue<block*> m_emptyWriteBuffers;
	std::queue<block*> m_emptyReadBuffers;
	std::queue<block*> m_fullWriteBuffers;
	std::queue<block*> m_fullReadBuffers;
	mutable boost::mutex m_mutex;
	mutable boost::condition_variable m_condition;

	leaves m_leaves;
	bits::tournament_tree<typename leaves::iterator, tourn_pred> m_tournamentTree;
};

} // namespace tpie

#endif // __TPIE_PIPELINING_MERGE_SORTER_H__
