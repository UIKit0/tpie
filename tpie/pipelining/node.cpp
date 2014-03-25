// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
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

#include <tpie/pipelining/node.h>

namespace tpie {

namespace pipelining {

/*virtual*/ node::~node() {}

memory_size_type node::get_minimum_memory() const {
	return m_minimumMemory;
}

memory_size_type node::get_maximum_memory() const {
	return m_maximumMemory;
}

memory_size_type node::get_available_memory() const {
	return m_availableMemory;
}

void node::set_memory_fraction(double f) {
	m_memoryFraction = f;
}

double node::get_memory_fraction() const {
	return m_memoryFraction;
}

bits::node_map::ptr node::get_node_map() const {
	return token.get_map();
}

node_token::id_t node::get_id() const {
	return token.id();
}

/*virtual*/ void node::prepare() {
}

/*virtual*/ void node::propagate() {
}

/*virtual*/ void node::begin() {
}

/*virtual*/ void node::go() {
	log_warning() << "node subclass " << typeid(*this).name()
		<< " is not an initiator node" << std::endl;
	throw not_initiator_node();
}

/*virtual*/ void node::end() {
}

/*virtual*/ bool node::can_evacuate() {
	return false;
}

/*virtual*/ void node::evacuate() {
}

priority_type node::get_name_priority() {
	return m_namePriority;
}

const std::string & node::get_name() {
	if (m_name.empty()) {
		m_name = bits::extract_pipe_name(typeid(*this).name());
	}
	return m_name;
}

void node::set_name(const std::string & name, priority_type priority /*= PRIORITY_USER*/) {
	m_name = name;
	m_namePriority = priority;
}

void node::set_breadcrumb(const std::string & breadcrumb) {
	m_name = m_name.empty() ? breadcrumb : (breadcrumb + " | " + m_name);
}

stream_size_type node::get_steps() {
	return m_stepsTotal;
}

void node::set_progress_indicator(progress_indicator_base * pi) {
	m_pi = pi;
}

progress_indicator_base * node::get_progress_indicator() {
	return m_pi;
}

node::STATE node::get_state() const {
	return m_state;
}

void node::set_state(STATE s) {
	m_state = s;
}

int node::get_plot_options() const {
	return m_plotOptions;
}

void node::set_plot_options(int options) {
	m_plotOptions = options;
}

node::node()
	: token(this)
	, m_minimumMemory(0)
	, m_maximumMemory(std::numeric_limits<memory_size_type>::max())
	, m_availableMemory(0)
	, m_memoryFraction(0.0)
	, m_namePriority(PRIORITY_NO_NAME)
	, m_stepsTotal(0)
	, m_stepsLeft(0)
	, m_pi(0)
	, m_state(STATE_FRESH)
	, m_plotOptions(0)
{
}

node::node(const node & other)
	: token(other.token, this)
	, m_minimumMemory(other.m_minimumMemory)
	, m_maximumMemory(other.m_maximumMemory)
	, m_availableMemory(other.m_availableMemory)
	, m_memoryFraction(other.m_memoryFraction)
	, m_name(other.m_name)
	, m_namePriority(other.m_namePriority)
	, m_stepsTotal(other.m_stepsTotal)
	, m_stepsLeft(other.m_stepsLeft)
	, m_pi(other.m_pi)
	, m_state(other.m_state)
	, m_plotOptions(other.m_plotOptions)
{
	if (m_state != STATE_FRESH) 
		throw call_order_exception(
			"Tried to copy pipeline node after prepare had been called");
}

#ifdef TPIE_CPP_RVALUE_REFERENCE
node::node(node && other)
	: token(std::move(other.token), this)
	, m_minimumMemory(std::move(other.m_minimumMemory))
	, m_maximumMemory(std::move(other.m_maximumMemory))
	, m_availableMemory(std::move(other.m_availableMemory))
	, m_memoryFraction(std::move(other.m_memoryFraction))
	, m_name(std::move(other.m_name))
	, m_namePriority(std::move(other.m_namePriority))
	, m_stepsTotal(std::move(other.m_stepsTotal))
	, m_stepsLeft(std::move(other.m_stepsLeft))
	, m_pi(std::move(other.m_pi))
	, m_state(std::move(other.m_state))
	, m_plotOptions(std::move(other.m_plotOptions))
{
	if (m_state != STATE_FRESH)
		throw call_order_exception(
			"Tried to move pipeline node after prepare had been called");
}
#endif // TPIE_CPP_RVALUE_REFERENCE

node::node(const node_token & token)
	: token(token, this, true)
	, m_minimumMemory(0)
	, m_maximumMemory(std::numeric_limits<memory_size_type>::max())
	, m_availableMemory(0)
	, m_memoryFraction(0.0)
	, m_namePriority(PRIORITY_NO_NAME)
	, m_stepsTotal(0)
	, m_stepsLeft(0)
	, m_pi(0)
	, m_state(STATE_FRESH)
	, m_plotOptions(0)
{
}

void node::add_push_destination(const node_token & dest) {
	bits::node_map::ptr m = token.map_union(dest);
	m->add_relation(token.id(), dest.id(), bits::pushes);
}

void node::add_push_destination(const node & dest) {
	if (get_state() != STATE_FRESH) {
		throw call_order_exception("add_push_destination called too late");
	}
	add_push_destination(dest.token);
}

void node::add_pull_source(const node_token & dest) {
	if (get_state() != STATE_FRESH) {
		throw call_order_exception("add_pull_source called too late");
	}
	bits::node_map::ptr m = token.map_union(dest);
	m->add_relation(token.id(), dest.id(), bits::pulls);
}

void node::add_pull_source(const node & dest) {
	add_pull_source(dest.token);
}

void node::add_dependency(const node_token & dest) {
	bits::node_map::ptr m = token.map_union(dest);
	m->add_relation(token.id(), dest.id(), bits::depends);
}

void node::add_dependency(const node & dest) {
	add_dependency(dest.token);
}

void node::set_minimum_memory(memory_size_type minimumMemory) {
	if (get_state() != STATE_FRESH && get_state() != STATE_IN_PREPARE) {
		throw call_order_exception("set_minimum_memory");
	}
	m_minimumMemory = minimumMemory;
}

void node::set_maximum_memory(memory_size_type maximumMemory) {
	if (get_state() != STATE_FRESH && get_state() != STATE_IN_PREPARE) {
		throw call_order_exception("set_maximum_memory");
	}
	m_maximumMemory = maximumMemory;
}

/*virtual*/ void node::set_available_memory(memory_size_type availableMemory) {
	m_availableMemory = availableMemory;
}

void node::forward_any(std::string key, boost::any value) {
	switch (get_state()) {
		case STATE_FRESH:
		case STATE_IN_PREPARE:
		case STATE_AFTER_PREPARE:
			// Allowed since forward() is allowed in prepare()
			break;
		case STATE_IN_PROPAGATE:
		case STATE_AFTER_PROPAGATE:
			// Allowed since forward() is allowed in propagate()
			break;
		case STATE_IN_BEGIN:
			throw call_order_exception("forward");
		case STATE_AFTER_BEGIN:
		case STATE_IN_END:
		case STATE_AFTER_END:
			// Allowed since forward() is allowed in end()
			break;
		default:
			log_debug() << "forward in unknown state " << get_state() << std::endl;
			break;
	}

	add_forwarded_data(key, value, true);

	bits::node_map::ptr nodeMap = get_node_map()->find_authority();

	typedef node_token::id_t id_t;
	std::vector<id_t> successors;
	nodeMap->get_successors(get_id(), successors);
	for (size_t i = 0; i < successors.size(); ++i) {
		nodeMap->get(successors[i])->add_forwarded_data(key, value, false);
	}
}

void node::add_forwarded_data(std::string key, boost::any value, bool explicitForward) {
	if (m_values.count(key) &&
		!explicitForward && m_values[key].second) return;
	m_values[key].first = value;
	m_values[key].second = explicitForward;
}

bool node::can_fetch(std::string key) {
	return m_values.count(key) != 0;
}

boost::any node::fetch_any(std::string key) {
	if (m_values.count(key) != 0) {
		return m_values[key].first;
	} else {
		std::stringstream ss;
		ss << "Tried to fetch nonexistent key '" << key << '\'';
		throw invalid_argument_exception(ss.str());
	}
}

const node_token & node::get_token() const {
	return token;
}

void node::set_steps(stream_size_type steps) {
	switch (get_state()) {
		case STATE_FRESH:
		case STATE_IN_PREPARE:
		case STATE_IN_PROPAGATE:
			break;
		case STATE_IN_BEGIN:
			log_error() << "set_steps in begin(); use set_steps in propagate() instead." << std::endl;
			throw call_order_exception("set_steps");
		default:
			log_error() << "set_steps in unknown state " << get_state() << std::endl;
			throw call_order_exception("set_steps");
	}
	m_stepsTotal = m_stepsLeft = steps;
}

void node::step(stream_size_type steps /*= 1*/) {
	assert(get_state() == STATE_IN_END || get_state() == STATE_AFTER_BEGIN || get_state() == STATE_IN_END);
	if (m_stepsLeft < steps) {
		if (m_stepsTotal != std::numeric_limits<stream_size_type>::max()) {
			log_warning() << typeid(*this).name() << " ==== Too many steps " << m_stepsTotal << std::endl;
			m_stepsLeft = 0;
			m_stepsTotal = std::numeric_limits<stream_size_type>::max();
		}
	} else {
		m_stepsLeft -= steps;
	}
	m_pi->step(steps);
}

progress_indicator_base * node::proxy_progress_indicator() {
	if (m_piProxy.get() != 0) return m_piProxy.get();
	progress_indicator_base * pi = new bits::proxy_progress_indicator(*this);
	m_piProxy.reset(pi);
	return pi;
}











} // namespace pipelining

} // namespace tpie
