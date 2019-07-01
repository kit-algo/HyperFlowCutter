#pragma once

#include "hla_flowcutter.h"
#include "../datastructure/queue.h"

namespace hyper {
	class HLAAssimilation {
	public:
		static void growAssimilated(const Hypergraph& hg, HLAFlowCutter& flc, LayeredQueue<nodeid>& queue) {
			queue.clear();
			for (nodeid s : flc.source_piercing_nodes) { queue.push(s); assert(!flc.isTarget(s)); assert(flc.isSource(s)); assert(!flc.target_nodes.contains(s)); }
			while (!queue.empty()) {
				nodeid u = queue.pop();
				assert(flc.isSource(u));
				assert(!flc.target_nodes.contains(u));
				for (netid e : hg.hyperedgesOf(u)) {
					if (!flc.source_hyperedges.is_assimilated(e)) {
						if (flc.flow_from[e] != u) {
							if (!flc.has_flow(e) || flc.flow_to[e] == u) {
								flc.hyperedge_assimilate(e);
								for (nodeid pin : hg.pinsOf(e)) {
									assert(!flc.target_nodes.contains(pin));
									if (!flc.isSource(pin)) {
										flc.node_assimilate(pin);
										queue.push(pin);
									}
								}
							}
							else {
								nodeid v = flc.flow_from[e];
								assert(flc.has_flow(e)); assert(flc.flow_from[e] != invalid_node);
								assert(!flc.target_nodes.contains(v));
								if (!flc.isSource(v)) {
									flc.node_assimilate(v);
									queue.push(v);
								}
								if (flc.shouldBeAddedToCutfront(e)) { flc.addToCut(e); }    //in theory the condition would not even need to be checked :). at the other branches all pins get assimilated :)
							}
						}
						else if (flc.shouldBeAddedToCutfront(e)) { flc.addToCut(e); }
					}
				}
			}
		}
	};


}