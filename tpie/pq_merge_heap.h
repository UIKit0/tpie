#ifndef _TPIE_PQ_MERGE_HEAP_H_
#define _TPIE_PQ_MERGE_HEAP_H_

#include "ami.h"
#include "tpie_log.h"
#include <cassert>

namespace tpie{

/////////////////////////////////////////////////////////
///
/// \class pq_merge_heap
/// \author Lars Hvam Petersen
///
/// pq_merge_heap
///
/////////////////////////////////////////////////////////
template<typename T, typename Comparator = std::less<T> >
class pq_merge_heap {
	public:
		/////////////////////////////////////////////////////////
		///
		/// Constructor
		///
		/// \param elements Maximum allowed size of the heap
		///
		/////////////////////////////////////////////////////////
		pq_merge_heap(TPIE_OS_SIZE_T elements);

		/////////////////////////////////////////////////////////
		///
		/// Destructor
		///
		/////////////////////////////////////////////////////////
		~pq_merge_heap();

		/////////////////////////////////////////////////////////
		///
		/// Insert an element into the priority queue
		///
		/// \param x The item
		/// \param run Where it comes from
		///
		/////////////////////////////////////////////////////////
		void push(const T& x, TPIE_OS_SIZE_T run);

		/////////////////////////////////////////////////////////
		///
		/// Remove the top element from the priority queue
		///
		/////////////////////////////////////////////////////////
		void pop();

		/////////////////////////////////////////////////////////
		///
		/// Remove the top element from the priority queue and 
		/// insert another
		///
		/// \param x The item
		/// \param run Where it comes from
		///
		/////////////////////////////////////////////////////////
		void pop_and_push(const T& x, TPIE_OS_SIZE_T run);

		/////////////////////////////////////////////////////////
		///
		/// See whats on the top of the priority queue
		///
		/// \return Top element
		///
		/////////////////////////////////////////////////////////
		const T& top();

		/////////////////////////////////////////////////////////
		///
		/// Return top element run number
		///
		/// \return Top element run number
		///
		/////////////////////////////////////////////////////////
		const TPIE_OS_SIZE_T top_run();

		/////////////////////////////////////////////////////////
		///
		/// Returns the size of the queue
		///
		/// \return Queue size
		///
		/////////////////////////////////////////////////////////
		const TPIE_OS_SIZE_T size();

		/////////////////////////////////////////////////////////
		///
		/// Return true if queue is empty otherwise false
		///
		/// \return Boolean - empty or not
		///
		/////////////////////////////////////////////////////////
		const bool empty();

	private:
		void fixDown();
		void validate();
		void dump();

		TPIE_OS_SIZE_T m_size;
		T min;
		Comparator comp_;

		T* heap;
		TPIE_OS_SIZE_T* runs;
		TPIE_OS_SIZE_T maxsize;
};

#include "pq_merge_heap.inl"
}
#endif
