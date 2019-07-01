#pragma once

#include "hla_flowcutter.h"
#include "../datastructure/queue.h"


namespace hyper{

class HLAGrowReachable {
public:
	static void growReachable(const Hypergraph& hg, HLAFlowCutter& flc, LayeredQueue<nodeid>& queue) {
		flc.source_nodes.reset_all_except_assimilated();
		flc.source_hyperedges.reset_all_except_assimilated();
		growReachableWithoutReset(hg, flc, queue);
	}
	static void growReachableWithoutReset(const Hypergraph& hg, HLAFlowCutter& flc, LayeredQueue<nodeid>& queue) {
		queue.clear();
		for (nodeid s : flc.source_piercing_nodes) { queue.push(s); assert(!flc.isTarget(s)); assert(flc.isSource(s)); assert(!flc.target_nodes.contains(s)); }
		while (!queue.empty()) {
			nodeid u = queue.pop();
			assert(flc.source_nodes.contains(u));
			assert(!flc.target_nodes.contains(u));
			for (netid e : hg.hyperedgesOf(u)) {
				if (!flc.source_hyperedges.contains(e) && flc.flow_from[e] != u) {
					if (!flc.has_flow(e) || flc.flow_to[e] == u) {
						flc.source_hyperedges.add_but_dont_assimilate(e);
						for (nodeid pin : hg.pinsOf(e)) {
							assert(!flc.target_nodes.contains(pin));
							if (!flc.source_nodes.contains(pin)) {
								flc.source_nodes.add_but_dont_assimilate(pin);
								queue.push(pin);
							}
						}
					}
					else {
						nodeid v = flc.flow_from[e];
						assert(!flc.target_nodes.contains(v));
						if (!flc.source_nodes.contains(v)) {
							flc.source_nodes.add_but_dont_assimilate(v);
							queue.push(v);
						}
					}
				}
			}
		}
	}
};
}