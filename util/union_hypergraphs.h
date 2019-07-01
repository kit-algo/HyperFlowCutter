#pragma once

#include "../datastructure/hypergraph.h"
namespace hyper {
class UnionDisjointHypergraphs {
public:
	template<typename Predicate>
	static Hypergraph unionDisjointHypergraphs(const std::vector<Hypergraph>& hgs, Predicate pred) {
		std::vector<nodeid> pins, hyperedge_sizes;
		nodeid global_node_id_offset = 0;
		auto global_id = [&global_node_id_offset](nodeid local) { return global_node_id_offset + local; };
		nodeid i = 0;
		for (auto& hg : hgs) {
			if (pred(i)) {
				for (netid e = 0; e < hg.numHyperedges(); e++) {
					hyperedge_sizes.push_back(hg.pinCount(e));
					std::transform(hg.pinsOf(e).begin(), hg.pinsOf(e).end(), std::back_inserter(pins), global_id);
				}
				global_node_id_offset += hg.totalNodeWeight();
			}
			i++;
		}
		return Hypergraph(global_node_id_offset, hyperedge_sizes, std::move(pins));
	}
};
}