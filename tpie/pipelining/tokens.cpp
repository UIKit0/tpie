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

#include <tpie/pipelining/tokens.h>
#include <tpie/pipelining/pipe_segment.h>
#include <tpie/pipelining/segment_base.h>

namespace tpie {

namespace pipelining {

segment_map::id_t segment_map::nextId = 0;

// Called by graph_traits
void segment_map::send_successors() const {
	for (relmapit i = m_relations.begin(); i != m_relations.end(); ++i) {
		if (i->second.second == uses) continue;
		pipe_segment * lhs = m_tokens.find(i->first)->second->assert_pipe_segment();
		pipe_segment * rhs = m_tokens.find(i->second.first)->second->assert_pipe_segment();
		switch (i->second.second) {
			case pushes:
				lhs->add_successor(rhs);
				break;
			case pulls:
			case depends:
				rhs->add_successor(lhs);
				break;
			case uses:
				break;
		}
	}
}

void segment_map::link(segment_map::ptr target) {
	if (target.get() == this) {
		// self link attempted
		// we must never have some_map->m_authority point to some_map,
		// since it would create a reference cycle
		return;
	}
	// union by rank
	if (target->m_rank > m_rank)
		return target->link(ptr(self));

	for (mapit i = target->begin(); i != target->end(); ++i) {
		set_token(i->first, i->second);
	}
	for (dsmapit i = target->m_dataStructures.begin(); i != target->m_dataStructures.end(); ++i) {
		set_data_structure(i->first, i->second);
	}
	for (relmapit i = target->m_relations.begin(); i != target->m_relations.end(); ++i) {
		m_relations.insert(*i);
	}
	for (relmapit i = target->m_relationsInv.begin(); i != target->m_relationsInv.end(); ++i) {
		m_relationsInv.insert(*i);
	}
	target->m_tokens.clear();
	target->m_authority = ptr(self);

	// union by rank
	if (target->m_rank == m_rank)
		++m_rank;
}

segment_map::ptr segment_map::find_authority() {
	if (!m_authority)
		return ptr(self);

	segment_map * i = m_authority.get();
	while (i->m_authority) {
		i = i->m_authority.get();
	}
	ptr result(i->self);

	// path compression
	segment_map * j = m_authority.get();
	while (j->m_authority) {
		segment_map * k = j->m_authority.get();
		j->m_authority = result;
		j = k;
	}

	return result;
}

size_t segment_map::out_degree(const relmap_t & map, id_t from, segment_relation rel) const {
	size_t res = 0;
	relmapit i = map.find(from);
	while (i != map.end() && i->first == from) {
		if (i->second.second == rel) ++res;
		++i;
	}
	return res;
}

segment_map::val_t segment_map::get(const segment_token & token) const {
	return get(token.id());
}

} // namespace pipelining

} // namespace tpie
