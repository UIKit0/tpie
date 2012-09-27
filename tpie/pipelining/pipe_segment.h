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

#ifndef __TPIE_PIPELINING_PIPE_SEGMENT_H__
#define __TPIE_PIPELINING_PIPE_SEGMENT_H__

#include <tpie/pipelining/exception.h>
#include <tpie/pipelining/tokens.h>
#include <tpie/progress_indicator_base.h>
#include <tpie/progress_indicator_null.h>
#include <boost/any.hpp>
#include <tpie/pipelining/data_structure.h>

namespace tpie {

namespace pipelining {

struct pipe_segment;

namespace bits {

class proxy_progress_indicator : public tpie::progress_indicator_base {
	pipe_segment & m_segment;

public:
	proxy_progress_indicator(pipe_segment & s)
		: progress_indicator_base(1)
		, m_segment(s)
	{
	}

	inline void refresh();
};

} // namespace bits

///////////////////////////////////////////////////////////////////////////////
/// Base class of all segments. A segment should inherit from pipe_segment,
/// have a single template parameter dest_t if it is not a terminus segment,
/// and implement methods begin(), push() and end(), if it is not a source
/// segment.
///////////////////////////////////////////////////////////////////////////////
struct pipe_segment : public segment_base {
	virtual void begin() {
		forward_all();
	}

	virtual void go() {
		progress_indicator_null pi;
		go(pi);
		// if go didn't throw, it was overridden - but it shouldn't be
		log_warning() << "pipe_segment subclass " << typeid(*this).name() << " uses old go() interface" << std::endl;
	}

	// Overriding this method is deprecated
	virtual void go(progress_indicator_base &) {
		log_warning() << "pipe_segment subclass " << typeid(*this).name() << " is not an initiator segment" << std::endl;
		throw not_initiator_segment();
	}

	virtual void end() {
	}

	virtual bool can_evacuate() {
		return false;
	}

	virtual void evacuate() {
	}

	// Called by segment_map
	inline void add_successor(pipe_segment * succ) {
		m_successors.push_back(succ);
	}

	inline stream_size_type get_steps() {
		return m_stepsTotal;
	}

	inline void set_progress_indicator(progress_indicator_base * pi) {
		m_pi = pi;
	}

protected:
	inline pipe_segment()
		: segment_base()
		, m_stepsTotal(0)
		, m_stepsLeft(0)
		, m_pi(0)
	{
		m_selfPipeSegment = this;
		m_selfDataStructure = 0;
	}

	inline pipe_segment(const pipe_segment & other)
		: segment_base(other)
		, m_stepsTotal(other.m_stepsTotal)
		, m_stepsLeft(other.m_stepsLeft)
		, m_pi(other.m_pi)
	{
		m_selfPipeSegment = this;
		m_selfDataStructure = 0;
	}

	inline pipe_segment(const segment_token & token)
		: segment_base(token)
		, m_stepsTotal(0)
		, m_stepsLeft(0)
		, m_pi(0)
	{
		m_selfPipeSegment = this;
		m_selfDataStructure = 0;
	}

	inline void add_push_destination(const segment_token & dest) {
		segment_map::ptr m = token.map_union(dest);
		m->add_relation(token.id(), dest.id(), pushes);
	}

	inline void add_push_destination(const pipe_segment & dest) {
		add_push_destination(dest.token);
	}

	inline void add_pull_destination(const segment_token & dest) {
		segment_map::ptr m = token.map_union(dest);
		m->add_relation(token.id(), dest.id(), pulls);
	}

	inline void add_pull_destination(const pipe_segment & dest) {
		add_pull_destination(dest.token);
	}

	inline void add_dependency(const segment_token & dest) {
		segment_map::ptr m = token.map_union(dest);
		m->add_relation(token.id(), dest.id(), depends);
	}

	inline void add_dependency(const pipe_segment & dest) {
		add_dependency(dest.token);
	}

	inline void add_data_structure(const data_structure & dest) {
		segment_map::ptr m = token.map_union(dest.token);
		m->add_relation(token.id(), dest.token.id(), uses);
	}

	template <typename T>
	inline void forward(std::string key, T value) {
		for (size_t i = 0; i < m_successors.size(); ++i) {
			m_successors[i]->m_values[key] = value;
		}
	}

	inline void forward_all() {
		for (valuemap::iterator i = m_values.begin(); i != m_values.end(); ++i) {
			forward(i->first, i->second);
		}
	}

	inline bool can_fetch(std::string key) {
		return m_values.count(key) != 0;
	}

	inline boost::any fetch_any(std::string key) {
		return m_values[key];
	}

	template <typename T>
	inline T fetch(std::string key) {
		return boost::any_cast<T>(m_values[key]);
	}

	const segment_token & get_token() {
		return token;
	}

	void set_steps(stream_size_type steps) {
		m_stepsTotal = m_stepsLeft = steps;
	}

	void step(stream_size_type steps = 1) {
		if (m_stepsLeft < steps) {
			log_warning() << typeid(*this).name() << " ==== Too many steps!" << std::endl;
			m_stepsLeft = 0;
		} else {
			m_stepsLeft -= steps;
		}
		m_pi->step(steps);
	}

	progress_indicator_base * proxy_progress_indicator() {
		if (m_piProxy.get() != 0) return m_piProxy.get();
		progress_indicator_base * pi = new bits::proxy_progress_indicator(*this);
		m_piProxy.reset(pi);
		return pi;
	}

	friend class phase;

private:
	std::vector<pipe_segment *> m_successors;
	typedef std::map<std::string, boost::any> valuemap;
	valuemap m_values;

	stream_size_type m_stepsTotal;
	stream_size_type m_stepsLeft;
	progress_indicator_base * m_pi;
	std::auto_ptr<progress_indicator_base> m_piProxy;

	friend class bits::proxy_progress_indicator;
};

namespace bits {

void proxy_progress_indicator::refresh() {
	double proxyMax = static_cast<double>(get_range());
	double proxyCur = static_cast<double>(get_current());
	double parentMax = static_cast<double>(m_segment.m_stepsTotal);
	double parentCur = static_cast<double>(m_segment.m_stepsTotal-m_segment.m_stepsLeft);
	double missing = parentMax*proxyCur/proxyMax - parentCur;
	if (missing < 1.0) return;
	stream_size_type times = static_cast<stream_size_type>(1.0+missing);
	times = std::min(m_segment.m_stepsLeft, times);
	m_segment.step(times);
}

} // namespace bits

} // namespace pipelining

} // namespace tpie

#endif // __TPIE_PIPELINING_PIPE_SEGMENT_H__
