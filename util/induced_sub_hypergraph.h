#pragma once
#include "../datastructure/hypergraph.h"

namespace hyper {
	class InducedSubHypergraph {
	public:
		template<typename VertexSubsetPredicate>
		static std::pair<Hypergraph, std::vector<nodeid>> extractInducedSubHypergraph(const Hypergraph& hg, VertexSubsetPredicate in_subset, bool keep_outside_hyperedges=false) {
			std::vector<nodeid> hyperedge_sizes, pins;
			std::vector<nodeid> __global2local(hg.numNodes(), invalid_node);
			std::vector<nodeid> __local2global;
			auto global2local = [&__global2local](nodeid x) { return __global2local[x]; };
			for (nodeid u = 0; u < hg.numNodes(); u++) {
				if (in_subset(u)) {
					__global2local[u] = static_cast<nodeid>(__local2global.size());
					__local2global.push_back(u);
				}
			}
			for (netid e = 0; e < hg.numHyperedges(); e++) {
				if (!keep_outside_hyperedges && !std::all_of(hg.pinsOf(e).begin(), hg.pinsOf(e).end(), in_subset)) {
					continue;
				}
				//filter, then transform
				hyperedge_sizes.push_back(0);
				for (nodeid pin : hg.pinsOf(e)) {
					if (in_subset(pin)) {
						hyperedge_sizes.back()++;
						pins.push_back(global2local(pin));
					}
				}
				if (hyperedge_sizes.back() == 0) { hyperedge_sizes.pop_back(); }
			}
			return std::make_pair(Hypergraph(static_cast<nodeid>(__local2global.size()), hyperedge_sizes, std::move(pins)), __local2global);
		}
	};
}