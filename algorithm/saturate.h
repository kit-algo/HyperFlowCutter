#pragma once


#include "hla_flowcutter.h"

namespace hyper {

	class Saturate {
	public:
		static inline void saturate(HLAFlowCutter& flc, nodeid u, netid e, nodeid v) {
			if (!flc.has_flow(e)) {
				//forward arc
				//flc.has_flow.set(e);
				flc.flow_from[e] = u; flc.flow_to[e] = v;
			}
			else if (flc.flow_from[e] == v && flc.flow_to[e] == u) {
				//back arc
				//flc.has_flow.reset(e);
				flc.flow_from[e] = invalid_node; flc.flow_to[e] = invalid_node;
			}
			else if (flc.flow_to[e] == u && flc.flow_from[e] != v) {
				flc.flow_to[e] = v;
			}
			else if (flc.flow_to[e] != u && flc.flow_from[e] == v) {
				flc.flow_from[e] = u;
			}
		}
	};




}