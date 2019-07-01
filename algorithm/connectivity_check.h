#pragma once

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "../datastructure/hypergraph.h"
#include <functional>

namespace hyper {

	//TODO rewrite to be generic
	class ConnectivityCheck {
	public:
		static bool is_connected(Hypergraph& hg) {
			std::vector<nodeid> queue(hg.numNodes());
			nodeid qfront = 0; nodeid qend = 1;
			queue[0] = 0;
			boost::dynamic_bitset<> was_node_seen(hg.numNodes());
			boost::dynamic_bitset<> was_hyperedge_seen(hg.numHyperedges());
			was_node_seen.set(0);
			nodeid num_nodes_seen = 1;
			while (qfront < qend && num_nodes_seen < hg.numNodes()) {
				nodeid u = queue[qfront++];
				for (netid e : hg.hyperedgesOf(u)) {
					if (!was_hyperedge_seen[e]) {
						was_hyperedge_seen.set(e);
						for (nodeid v : hg.pinsOf(e)) {
							if (!was_node_seen[v]) {
								num_nodes_seen++;
								was_node_seen.set(v);
								queue[qend++] = v;
							}
						}
					}
				}
			}
			return num_nodes_seen == hg.numNodes();
		}

		static std::vector<std::pair<uint32_t, uint32_t>> componentSizes(Hypergraph& hg) {
			if (hg.hasCoarseNodes()) throw std::runtime_error("Component sizes analytics not supported for coarse hypergraphs. Sorry.");
			std::vector<nodeid> ccsizes;

			std::vector<nodeid> queue(hg.numNodes());
			boost::dynamic_bitset<> was_node_seen(hg.numNodes());
			boost::dynamic_bitset<> was_hyperedge_seen(hg.numHyperedges());

			for (nodeid s = 0; s < hg.numNodes(); s++) {
				if (!was_node_seen[s]) {
					nodeid qfront = 0; nodeid qend = 1;
					queue[0] = s;
					was_node_seen.set(s);
					nodeid num_nodes_seen = 1;
					while (qfront < qend) {
						nodeid u = queue[qfront++];
						for (netid e : hg.hyperedgesOf(u)) {
							if (!was_hyperedge_seen[e]) {
								was_hyperedge_seen.set(e);
								for (nodeid v : hg.pinsOf(e)) {
									if (!was_node_seen[v]) {
										num_nodes_seen++;
										was_node_seen.set(v);
										queue[qend++] = v;
									}
								}
							}
						}
					}
					ccsizes.push_back(num_nodes_seen);
				}
			}
			std::sort(ccsizes.begin(), ccsizes.end(), std::greater<>());
			std::vector<std::pair<uint32_t, uint32_t>> histogram;
			for (auto& x : ccsizes) {
				if (histogram.empty() || histogram.back().first != x) {
					histogram.emplace_back(x, 1);
				}
				else {
					 histogram.back().second++;
				}
			}
			return histogram;
		}

		template<typename FScanCondition, typename FEnqueueCondition>
		static void BFSWithScanCondition(Hypergraph& hg, nodeid s, FScanCondition scan_condition, FEnqueueCondition qc) {
			std::vector<nodeid> queue(hg.numNodes());
			nodeid qfront = 0; nodeid qend = 1;
			queue[0] = s;
			boost::dynamic_bitset<> was_node_seen(hg.numNodes());
			boost::dynamic_bitset<> was_hyperedge_seen(hg.numHyperedges());
			was_node_seen.set(s);
			while (qfront < qend) {
				nodeid u = queue[qfront++];
				if (!scan_condition(u)) { continue; }
				for (netid e : hg.hyperedgesOf(u)) {
					if (!was_hyperedge_seen[e]) {
						was_hyperedge_seen.set(e);
						for (nodeid v : hg.pinsOf(e)) {
							if (!was_node_seen[v] && qc(v)) {
								was_node_seen.set(v);
								queue[qend++] = v;
							}
						}
					}
				}
			}
			std::cout << "Reached " << was_node_seen.count() << " vertices out of " << hg.numNodes() << std::endl;
		}

		static ConnectedComponents connectedComponents(Hypergraph& hg) {
			if (hg.hasCoarseNodes()) throw std::runtime_error("CCs not supported for coarse hypergraphs. Sorry.");
			std::vector<nodeid> node_component(hg.numNodes(), INF);
			std::vector<nodeid> hyperedge_component(hg.numHyperedges(), INF);
			std::vector<nodeid> ccsizes;

			std::vector<nodeid> queue(hg.numNodes());
			nodeid component_id = 0;
			for (nodeid s = 0; s < hg.numNodes(); s++) {
				if (node_component[s] == INF) {
					nodeid qfront = 0; nodeid qend = 1;
					queue[0] = s;
					node_component[s] = component_id;
					nodeid num_nodes_seen = 1;
					while (qfront < qend) {
						nodeid u = queue[qfront++];
						for (netid e : hg.hyperedgesOf(u)) {
							if (hyperedge_component[e] == INF) {
								hyperedge_component[e] = component_id;
								for (nodeid v : hg.pinsOf(e)) {
									if (node_component[v] == INF) {
										num_nodes_seen++;
										node_component[v] = component_id;
										queue[qend++] = v;
									}
								}
							}
						}
					}
					ccsizes.push_back(num_nodes_seen);
					component_id++;
				}
			}
			return { node_component, hyperedge_component, ccsizes };
		}
	};

}