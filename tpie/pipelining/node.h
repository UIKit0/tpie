// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :
// Copyright 2012, The TPIE development team
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

#ifndef __TPIE_PIPELINING_NODE_H__
#define __TPIE_PIPELINING_NODE_H__

#include <tpie/pipelining/exception.h>
#include <tpie/pipelining/tokens.h>
#include <tpie/progress_indicator_base.h>
#include <tpie/progress_indicator_null.h>
#include <boost/any.hpp>
#include <tpie/pipelining/priority_type.h>
#include <tpie/pipelining/predeclare.h>
#include <tpie/pipelining/node_name.h>
#ifdef TPIE_CPP_DECLTYPE
#include <type_traits>
#endif
#include <utility>

namespace tpie {

namespace pipelining {

namespace bits {

class proxy_progress_indicator : public tpie::progress_indicator_base {
	node & m_node;

public:
	proxy_progress_indicator(node & s)
		: progress_indicator_base(1)
		, m_node(s)
	{
	}

	inline void refresh();
};


template <typename T>
struct remove {
	typedef T type;
};

template <typename T>
struct remove<const T> {
	typedef typename remove<T>::type type;
};

template <typename T>
struct remove<T &> {
	typedef typename remove<T>::type type;
};

#ifdef TPIE_CPP_RVALUE_REFERENCE
template <typename T>
struct remove<T &&> {
	typedef typename remove<T>::type type;
};
#endif // TPIE_CPP_RVALUE_REFERENCE


#ifdef TPIE_CPP_DECLTYPE

template <typename T>
struct push_traits {};

template <typename ClassType, typename ArgType>
struct push_traits< void(ClassType::*)(ArgType) > {
	typedef ArgType type;
};

template <typename T, bool>
struct push_type_help {
	typedef typename push_traits<decltype(&T::push)>::type type_;
	typedef typename remove<type_>::type type;
};

template <typename T>
struct push_type_help<T, true> {
	typedef typename T::item_type type;
};

template <typename T>
struct pull_traits {};

template <typename ClassType, typename ArgType>
struct pull_traits< ArgType(ClassType::*)() > {
	typedef ArgType type;
};

template <typename T, bool>
struct pull_type_help {
	typedef typename pull_traits<decltype(&T::pull)>::type type_;
	typedef typename remove<type_>::type type;
};

template <typename T>
struct pull_type_help<T, true> {
	typedef typename T::item_type type;
};

template <typename T>
struct has_itemtype {
	typedef char yes[1];
	typedef char no[2];
	
	template <typename C>
	static yes& test(typename C::item_type*);
	
	template <typename>
	static no& test(...);
	//static_assert(sizeof(test<T>(nullptr)) == sizeof(yes), "WTF");
	static const bool value = sizeof(test<T>(nullptr)) == sizeof(yes);
};
#endif //TPIE_CPP_DECLTYPE


} // namespace bits

///////////////////////////////////////////////////////////////////////////////
/// Class to deduce the item_type of a node of type T.
/// typedef item_type to T::item_type if that exists
/// otherwise typedef it to the type of whatever paramater push takes
///////////////////////////////////////////////////////////////////////////////
#ifdef TPIE_CPP_DECLTYPE
template <typename T>
struct push_type {
	typedef typename bits::remove<T>::type node_type;
	typedef typename bits::push_type_help<node_type, bits::has_itemtype<T>::value>::type type;
};

template <typename T>
struct pull_type {
	typedef typename bits::remove<T>::type node_type;
	typedef typename bits::pull_type_help<node_type, bits::has_itemtype<T>::value>::type type;
};

#else //TPIE_CPP_DECLTYPE

template <typename T>
struct push_type {
  typedef typename T::item_type type;
};

template <typename T>
struct pull_type {
  typedef typename T::item_type type;
};

#endif //TPIE_CPP_DECLTYPE

///////////////////////////////////////////////////////////////////////////////
/// Base class of all nodes. A node should inherit from the node class,
/// have a single template parameter dest_t if it is not a terminus node,
/// and implement methods begin(), push() and end(), if it is not a source
/// node.
///////////////////////////////////////////////////////////////////////////////
class node {
public:
	///////////////////////////////////////////////////////////////////////////
	// \brief Options for how to plot this node
	//////////////////////////////////////////////////////////////////////////
	enum PLOT {
		PLOT_SIMPLIFIED_HIDE=1,
		PLOT_BUFFERED=2
	};

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Used internally to check order of method calls.
	///////////////////////////////////////////////////////////////////////////
	enum STATE {
		STATE_FRESH,
		STATE_IN_PREPARE,
		STATE_AFTER_PREPARE,
		STATE_IN_PROPAGATE,
		STATE_AFTER_PROPAGATE,
		STATE_IN_BEGIN,
		STATE_AFTER_BEGIN,
		STATE_IN_END,
		STATE_AFTER_END
	};

	///////////////////////////////////////////////////////////////////////////
	/// \brief Virtual dtor.
	///////////////////////////////////////////////////////////////////////////
	virtual ~node();

	///////////////////////////////////////////////////////////////////////////
	/// \brief Get the minimum amount of memory declared by this node.
	/// Defaults to zero when no minimum has been set.
	///////////////////////////////////////////////////////////////////////////
	memory_size_type get_minimum_memory() const;

	///////////////////////////////////////////////////////////////////////////
	/// \brief Get the maximum amount of memory declared by this node.
	/// Defaults to maxint when no maximum has been set.
	///////////////////////////////////////////////////////////////////////////
	memory_size_type get_maximum_memory() const;

	///////////////////////////////////////////////////////////////////////////
	/// \brief Get the amount of memory assigned to this node.
	///////////////////////////////////////////////////////////////////////////
	memory_size_type get_available_memory() const;

	///////////////////////////////////////////////////////////////////////////
	/// \brief Set the memory priority of this node. Memory is distributed
	/// proportionally to the priorities of the nodes in the given phase.
	///////////////////////////////////////////////////////////////////////////
	void set_memory_fraction(double f);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Get the memory priority of this node.
	///////////////////////////////////////////////////////////////////////////
	double get_memory_fraction() const;

	///////////////////////////////////////////////////////////////////////////
	/// \brief Get the local node map, mapping node IDs to node
	/// pointers for all the nodes reachable from this one.
	///////////////////////////////////////////////////////////////////////////
	bits::node_map::ptr get_node_map() const;

	///////////////////////////////////////////////////////////////////////////
	/// \brief Get the internal node ID of this node (mainly
	/// for debugging purposes).
	///////////////////////////////////////////////////////////////////////////
	node_token::id_t get_id() const;

	///////////////////////////////////////////////////////////////////////////
	/// \brief Called before memory assignment but after depending phases have
	/// executed and ended. The implementer may use fetch and forward in this
	/// phase. The implementer does not have to call the super prepare-method;
	/// its default implementation is empty.
	///////////////////////////////////////////////////////////////////////////
	virtual void prepare();

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Propagate stream metadata.
	///
	/// The implementation may fetch() and forward() metadata such as number of
	/// items or the size of a single item.
	///
	/// The pipelining framework calls propagate() on the nodes in the
	/// item flow graph in a topological order.
	///
	/// The default implementation does nothing.
	///////////////////////////////////////////////////////////////////////////
	virtual void propagate();

	///////////////////////////////////////////////////////////////////////////
	/// \brief Begin pipeline processing phase.
	///
	/// The implementation may pull() from a pull destination in begin(),
	/// and it may push() to a push destination.
	///
	/// The pipelining framework calls begin() on the nodes in the
	/// actor graph in a reverse topological order. The framework calls
	/// node::begin() on a node after calling begin() on its pull and push
	/// destinations.
	///
	/// The default implementation does nothing.
	///////////////////////////////////////////////////////////////////////////
	virtual void begin();

	///////////////////////////////////////////////////////////////////////////
	/// \brief For initiator nodes, execute this phase by pushing all items
	/// to be pushed. For non-initiator nodes, the default implementation
	/// throws a not_initiator_segment exception.
	///////////////////////////////////////////////////////////////////////////
	virtual void go();

	///////////////////////////////////////////////////////////////////////////
	/// \brief End pipeline processing phase.
	///
	/// The implementation may pull() from a pull destination in end(),
	/// and it may push() to a push destination.
	///
	/// The pipelining framework calls end() on the nodes in the
	/// pipeline graph in a topological order. The framework calls
	/// node::end() on a node before its pull and push
	/// destinations.
	///
	/// The default implementation does nothing, so it does not matter if the
	/// implementation calls the parent end().
	///////////////////////////////////////////////////////////////////////////
	virtual void end();

	///////////////////////////////////////////////////////////////////////////
	/// \brief Overridden by nodes that have data to evacuate.
	///////////////////////////////////////////////////////////////////////////
	virtual bool can_evacuate();

	///////////////////////////////////////////////////////////////////////////
	/// \brief Overridden by nodes that have data to evacuate.
	///////////////////////////////////////////////////////////////////////////
	virtual void evacuate();

	///////////////////////////////////////////////////////////////////////////
	/// \brief Get the priority of this node's name. For purposes of
	/// pipeline debugging and phase naming for progress indicator breadcrumbs.
	///////////////////////////////////////////////////////////////////////////
	priority_type get_name_priority();

	///////////////////////////////////////////////////////////////////////////
	/// \brief Get this node's name. For purposes of pipeline debugging and
	/// phase naming for progress indicator breadcrumbs.
	///////////////////////////////////////////////////////////////////////////
	const std::string & get_name();

	///////////////////////////////////////////////////////////////////////////
	/// \brief Set this node's name. For purposes of pipeline debugging and
	/// phase naming for progress indicator breadcrumbs.
	///////////////////////////////////////////////////////////////////////////
	void set_name(const std::string & name, priority_type priority = PRIORITY_USER);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Used internally when a pair_factory has a name set.
	///////////////////////////////////////////////////////////////////////////
	void set_breadcrumb(const std::string & breadcrumb);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Used internally for progress indication. Get the number of times
	/// the node expects to call step() at most.
	///////////////////////////////////////////////////////////////////////////
	stream_size_type get_steps();

	///////////////////////////////////////////////////////////////////////////
	/// \brief Used internally. Set the progress indicator to use.
	///////////////////////////////////////////////////////////////////////////
	void set_progress_indicator(progress_indicator_base * pi);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Used internally. Get the progress indicator used.
	///////////////////////////////////////////////////////////////////////////
	progress_indicator_base * get_progress_indicator();

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Used internally to check order of method calls.
	///////////////////////////////////////////////////////////////////////////
	STATE get_state() const;

	///////////////////////////////////////////////////////////////////////////
	/// \brief  Used internally to check order of method calls.
	///////////////////////////////////////////////////////////////////////////
	void set_state(STATE s);

	int get_plot_options() const;

	void set_plot_options(int options);
protected:
#ifdef _WIN32
	// Disable warning C4355: 'this' : used in base member initializer list
	// node_token does not access members of the `node *`,
	// it merely uses it as a value in the node map.
	// Only after this node object is completely constructed are node members accessed.
#pragma warning( push )
#pragma warning( disable : 4355 )
#endif // _WIN32
	///////////////////////////////////////////////////////////////////////////
	/// \brief Default constructor, using a new node_token.
	///////////////////////////////////////////////////////////////////////////
	node();

	///////////////////////////////////////////////////////////////////////////
	/// \brief Copy constructor. We need to define this explicitly since the
	/// node_token needs to know its new owner.
	///////////////////////////////////////////////////////////////////////////
	node(const node & other);

#ifdef TPIE_CPP_RVALUE_REFERENCE
	///////////////////////////////////////////////////////////////////////////
	/// \brief Move constructor. We need to define this explicitly since the
	/// node_token needs to know its new owner.
	///////////////////////////////////////////////////////////////////////////
	node(node && other);
#endif // TPIE_CPP_RVALUE_REFERENCE

	///////////////////////////////////////////////////////////////////////////
	/// \brief Constructor using a given fresh node_token.
	///////////////////////////////////////////////////////////////////////////
	node(const node_token & token);
#ifdef _WIN32
#pragma warning( pop )
#endif // _WIN32

	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by implementers to declare a push destination.
	///////////////////////////////////////////////////////////////////////////
	void add_push_destination(const node_token & dest);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by implementers to declare a push destination.
	///////////////////////////////////////////////////////////////////////////
	void add_push_destination(const node & dest);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by implementers to declare a pull source.
	///////////////////////////////////////////////////////////////////////////
	void add_pull_source(const node_token & dest);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by implementers to declare a pull source.
	///////////////////////////////////////////////////////////////////////////
	void add_pull_source(const node & dest);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by implementers to declare a node dependency, that is,
	/// a requirement that another node has end() called before the begin()
	/// of this node.
	///////////////////////////////////////////////////////////////////////////
	void add_dependency(const node_token & dest);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by implementers to declare a node dependency, that is,
	/// a requirement that another node has end() called before the begin()
	/// of this node.
	///////////////////////////////////////////////////////////////////////////
	void add_dependency(const node & dest);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by implementers to declare minimum memory requirements.
	///////////////////////////////////////////////////////////////////////////
	void set_minimum_memory(memory_size_type minimumMemory);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by implementers to declare maximum memory requirements.
	///
	/// To signal that you don't want any memory, set minimum memory and the
	/// memory fraction to zero.
	///////////////////////////////////////////////////////////////////////////
	void set_maximum_memory(memory_size_type maximumMemory);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by the memory manager to set the amount of memory
	/// assigned to this node.
	///////////////////////////////////////////////////////////////////////////
	virtual void set_available_memory(memory_size_type availableMemory);

public:
	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by implementers to forward auxiliary data to successors.
	/// If explicitForward is false, the data will not override data forwarded
	/// with explicitForward == true.
	///////////////////////////////////////////////////////////////////////////
	// Implementation note: If the type of the `value` parameter is changed
	// from `T` to `const T &`, this will yield linker errors if an application
	// attempts to pass a const reference to a static data member inside a
	// templated class.
	// See http://stackoverflow.com/a/5392050
	///////////////////////////////////////////////////////////////////////////
	template <typename T>
	void forward(std::string key, T value) {
		forward_any(key, boost::any(value));
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief See \ref node::forward.
	///////////////////////////////////////////////////////////////////////////
	void forward_any(std::string key, boost::any value);

private:
	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by forward_any to add forwarded data.
	//
	/// If explicitForward is false, the data will not override data forwarded
	/// with explicitForward == true.
	///////////////////////////////////////////////////////////////////////////
	void add_forwarded_data(std::string key, boost::any value, bool explicitForward);

public:
	///////////////////////////////////////////////////////////////////////////
	/// \brief Find out if there is a piece of auxiliary data forwarded with a
	/// given name.
	///////////////////////////////////////////////////////////////////////////
	bool can_fetch(std::string key);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Fetch piece of auxiliary data as boost::any (the internal
	/// representation).
	///////////////////////////////////////////////////////////////////////////
	boost::any fetch_any(std::string key);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Fetch piece of auxiliary data, expecting a given value type.
	///////////////////////////////////////////////////////////////////////////
	template <typename T>
	inline T fetch(std::string key) {
		if (m_values.count(key) == 0) {
			std::stringstream ss;
			ss << "Tried to fetch nonexistent key '" << key
			   << "' of type " << typeid(T).name();
			throw invalid_argument_exception(ss.str());
		}
		try {
			return boost::any_cast<T>(m_values[key].first);
		} catch (boost::bad_any_cast m) {
			std::stringstream ss;
			ss << "Trying to fetch key '" << key << "' of type "
			   << typeid(T).name() << " but forwarded data was of type "
			   << m_values[key].first.type().name() << ". Message was: " << m.what();
			throw invalid_argument_exception(ss.str());
		}
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Get the node_token that maps this node's ID to a pointer
	/// to this.
	///////////////////////////////////////////////////////////////////////////
	const node_token & get_token() const;

public:
	///////////////////////////////////////////////////////////////////////////
	/// \brief Called by implementers that intend to call step().
	/// \param steps  The number of times step() is called at most.
	///////////////////////////////////////////////////////////////////////////
	void set_steps(stream_size_type steps);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Step the progress indicator.
	/// \param steps  How many steps to step.
	///////////////////////////////////////////////////////////////////////////
	void step(stream_size_type steps = 1);

	///////////////////////////////////////////////////////////////////////////
	/// \brief Get a non-initialized progress indicator for use with external
	/// implementations. When step() is called on a proxying progress
	/// indicator, step() is called on the node according to the number
	/// of steps declared in progress_indicator_base::init() and in
	/// node::set_steps().
	///////////////////////////////////////////////////////////////////////////
	progress_indicator_base * proxy_progress_indicator();

#ifdef DOXYGEN
	///////////////////////////////////////////////////////////////////////////
	/// \brief For pull nodes, return true if there are more items to be
	/// pulled.
	///////////////////////////////////////////////////////////////////////////
	inline bool can_pull() const;

	///////////////////////////////////////////////////////////////////////////
	/// \brief For pull nodes, pull the next item from this node.
	///////////////////////////////////////////////////////////////////////////
	inline item_type pull();

	///////////////////////////////////////////////////////////////////////////
	/// \brief For push nodes, push the next item to this node.
	///////////////////////////////////////////////////////////////////////////
	inline void push(const item_type & item);
#endif

	friend class bits::memory_runtime;

	friend class factory_base;

private:
	node_token token;

	memory_size_type m_minimumMemory;
	memory_size_type m_maximumMemory;
	memory_size_type m_availableMemory;
	double m_memoryFraction;

	std::string m_name;
	priority_type m_namePriority;

	typedef std::map<std::string, std::pair<boost::any, bool> > valuemap;
	valuemap m_values;

	stream_size_type m_stepsTotal;
	stream_size_type m_stepsLeft;
	progress_indicator_base * m_pi;
	STATE m_state;
	std::auto_ptr<progress_indicator_base> m_piProxy;
	int m_plotOptions;

	friend class bits::proxy_progress_indicator;
};

namespace bits {

void proxy_progress_indicator::refresh() {
	double proxyMax = static_cast<double>(get_range());
	double proxyCur = static_cast<double>(get_current());
	double parentMax = static_cast<double>(m_node.m_stepsTotal);
	double parentCur = static_cast<double>(m_node.m_stepsTotal-m_node.m_stepsLeft);
	double missing = parentMax*proxyCur/proxyMax - parentCur;
	if (missing < 1.0) return;
	stream_size_type times = static_cast<stream_size_type>(1.0+missing);
	times = std::min(m_node.m_stepsLeft, times);
	m_node.step(times);
}

} // namespace bits

} // namespace pipelining

} // namespace tpie

#endif // __TPIE_PIPELINING_NODE_H__
