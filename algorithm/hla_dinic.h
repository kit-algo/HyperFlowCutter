#pragma once

#include "hypergraph_search_algorithms.h"
#include "hla_flowcutter.h"
#include "../datastructure/stack.h"
#include "../datastructure/queue.h"
#include "saturate.h"
#include "grow_assimilated.h"
#include "grow_reachable.h"

namespace hyper {
	class HLADinic : public SearchAlgorithm {
	private:
		static constexpr bool debug = false;
		static constexpr bool gather_stats = false;

		struct StackElement {
			nodeid tail;
		};

		using distance_t = nodeid;
		using pin_iterator = Hypergraph::pin_iterator;
		using hyperedge_iterator = Hypergraph::hyperedge_iterator;

		static constexpr distance_t inf_dist = invalid_node;
		static constexpr distance_t default_stale_distance = 0;
		static constexpr distance_t initial_base_distance = 2;

		FixedCapacityStack<nodeid, StackElement> stack;
		LayeredQueue<nodeid> queue;
		std::vector<pin_iterator> current_pin;
		std::vector<hyperedge_iterator> current_hyperedge;
		std::vector<distance_t> node_distance;
		std::vector<distance_t> hyperedge_distance;
		distance_t base_distance = initial_base_distance;
		distance_t bfs_dist = inf_dist;
		bool build_layered_network_during_grow_reachable;

	public:


		class CutterSpecificDatastructures : public SearchAlgorithm::CutterSpecificDatastructures {
		public:
			std::vector<pin_iterator> csd_current_pin;
			std::vector<hyperedge_iterator> csd_current_hyperedge;
			std::vector<distance_t> csd_node_distance;
			std::vector<distance_t> csd_hyperedge_distance;
			std::array<distance_t, 2> base_distance_next_bfs = { initial_base_distance, initial_base_distance };
			std::array<distance_t, 2> base_distance_previous_bfs = { initial_base_distance, initial_base_distance };

			std::array<bool, 2> layered_network_was_built = { false, false };

			explicit CutterSpecificDatastructures(const Hypergraph &hg) : csd_current_pin(hg.numHyperedges()),
																		  csd_current_hyperedge(hg.numNodes()),
																		  csd_node_distance(hg.numNodes(), default_stale_distance),
																		  csd_hyperedge_distance(hg.numHyperedges(), default_stale_distance)
			{
				for (nodeid u = 0; u < hg.numNodes(); u++)
					csd_current_hyperedge[u] = hg.endHyperedges(u);
				for (netid e = 0; e < hg.numHyperedges(); e++)
					csd_current_pin[e] = hg.endPins(e);
			}

			virtual CutterSpecificDatastructures* clone_impl() const override { return new CutterSpecificDatastructures(*this); }
		};


		explicit HLADinic(const Hypergraph& hg, bool build_datastructures_during_grow_reachable) :
				SearchAlgorithm(hg, build_datastructures_during_grow_reachable),
				stack(hg.numNodes()),
				queue(hg.numNodes()),
				build_layered_network_during_grow_reachable(build_datastructures_during_grow_reachable)
		{ }

		std::unique_ptr<SearchAlgorithm::CutterSpecificDatastructures> obtainOwnedCutterSpecificDatastructures() const override {
			return std::make_unique<CutterSpecificDatastructures>(hg);
		}

		void growAssimilated(HLAFlowCutter& flc) override { return HLAAssimilation::growAssimilated(hg, flc, queue); }

		void growReachableWithoutReset(HLAFlowCutter& flc) override {
			bool source_piercing_node_is_contracted_and_high_degree = false;
			if (flc.source_piercing_nodes.size() == 1) {
				nodeid s = flc.source_piercing_nodes[0];
				source_piercing_node_is_contracted_and_high_degree = (hg.nodeWeight(s) != 1) && (hg.degree(s) > 0.11 * hg.numHyperedges());
			}

			bool skip_building_layered_network_during_grow_reachable = source_piercing_node_is_contracted_and_high_degree; // || source_piercing_nodes_have_high_degree_sum;

			if (!build_layered_network_during_grow_reachable || skip_building_layered_network_during_grow_reachable)
				growReachableWithoutReset_Standard(flc);
			else
				growReachableWithoutReset_BuildingLayeredNetwork(flc);

		}

		void growReachableWithoutReset_Standard(HLAFlowCutter& flc) {
			HLAGrowReachable::growReachableWithoutReset(hg, flc, queue);
			static_cast<CutterSpecificDatastructures&>(*flc.csd).layered_network_was_built[flc.source_side_in_current_direction] = false;
		}

		void growReachableWithoutReset_BuildingLayeredNetwork(HLAFlowCutter& flc) {
			auto& csd = static_cast<CutterSpecificDatastructures&>(*flc.csd);
			swapCutterSpecificDatastructuresAndLocalDatastructures(csd);
			verify_datastructures_are_acquired();
			csd.base_distance_previous_bfs[flc.source_side_in_current_direction] = csd.base_distance_next_bfs[flc.source_side_in_current_direction];
			base_distance = csd.base_distance_next_bfs[flc.source_side_in_current_direction];
			growReachableAndBuildLayeredNetworkImpl(hg, flc);
			csd.base_distance_next_bfs[flc.source_side_in_current_direction] = bfs_dist;

			swapCutterSpecificDatastructuresAndLocalDatastructures(csd);
			verify_datastructures_are_released();
			csd.layered_network_was_built[flc.source_side_in_current_direction] = true;
		}

		void growReachable(HLAFlowCutter& flc) override {
			flc.source_nodes.reset_all_except_assimilated();
			flc.source_hyperedges.reset_all_except_assimilated();
			growReachableWithoutReset(flc);
		}

		flow_t exhaustFlow(HLAFlowCutter& flc) override {
			assert(flc.augmentingPathExists);
			auto& csd = static_cast<CutterSpecificDatastructures&>(*flc.csd);
			swapCutterSpecificDatastructuresAndLocalDatastructures(csd);
			verify_datastructures_are_acquired();
			flow_t flow = 0;
			flow_t recycle_flow = 0;
			flow_t pure_dinic_flow = 0;

			if (!flc.firstGrowReachable && csd.layered_network_was_built[1 - flc.source_side_in_current_direction] && build_layered_network_during_grow_reachable) {
				//flip. use dinicDFS. flip.
				auto t = time_now();
				flc.flipSearchDirection();
				base_distance = csd.base_distance_previous_bfs[flc.source_side_in_current_direction];
				recycle_flow += dinicDFS(hg, flc);
				flc.flipSearchDirection();
				LOG << V(recycle_flow) << inMilliseconds(duration(time_now() -t )).count() << "ms for recycle flow";
			}

			base_distance = std::max(
					csd.base_distance_next_bfs[flc.source_side_in_current_direction],
					csd.base_distance_next_bfs[1 - flc.source_side_in_current_direction]
			);

			auto t = time_now();
			while (dinicBFS(hg, flc)) {
				flow_t flow_diff = dinicDFS(hg, flc);
				pure_dinic_flow += flow_diff;
				base_distance = bfs_dist;
				LOG << V(flow_diff) << V(pure_dinic_flow);
			}

			flow = recycle_flow + pure_dinic_flow;
			LOG << V(flow) << inMilliseconds(duration(time_now() - t)).count() << "ms for regular dinic";

			verify_flow_is_exhausted(hg, flc);

			{
				csd.layered_network_was_built[flc.source_side_in_current_direction] = true;
				csd.base_distance_previous_bfs[flc.source_side_in_current_direction] = base_distance;
				csd.base_distance_next_bfs[flc.source_side_in_current_direction] = bfs_dist;
				csd.base_distance_next_bfs[1 - flc.source_side_in_current_direction] = base_distance;
				csd.base_distance_previous_bfs[1-flc.source_side_in_current_direction] = base_distance;
				swapCutterSpecificDatastructuresAndLocalDatastructures(csd);
				verify_datastructures_are_released();
			}
			return flow;
		}

	private:
		inline bool is_base_distance_safe() const {
			return static_cast<std::size_t>(base_distance) + static_cast<std::size_t>(hg.numNodes()) + 5 <= static_cast<std::size_t>(std::numeric_limits<distance_t>::max());
		}
		inline bool is_node_distance_stale(nodeid u) const { return node_distance[u] < base_distance; }
		inline bool is_hyperedge_distance_stale(netid e) const { return hyperedge_distance[e] < base_distance; }
		inline bool reset_distances_if_necessary() {
			verify_distances_are_stale();
			if (!is_base_distance_safe()) {
				LOG << "Resetting distance array" << V(base_distance);
				node_distance.assign(node_distance.size(), default_stale_distance);
				hyperedge_distance.assign(hyperedge_distance.size(), default_stale_distance);
				base_distance = initial_base_distance;
				return true;
			}
			return false;
		}

		template<typename FDontVisitAndMarkNodeIf, typename FQueueStoppingCriterion>
		void buildLayeredNetwork(const Hypergraph& hg, HLAFlowCutter& flc, FDontVisitAndMarkNodeIf dont_visit_and_mark_node_if, FQueueStoppingCriterion stop_search) {
			assert(is_base_distance_safe());
			bfs_dist = base_distance;
			queue.clear();

			auto visit_node = [&](nodeid v) {
				if (gather_stats) { bfsstats.visited_nodes++; }
				assert(!flc.isTarget(v));
				assert(is_node_distance_stale(v));
				node_distance[v] = bfs_dist;
				current_hyperedge[v] = hg.beginHyperedges(v);
				queue.push(v);
			};

			auto visit_node_and_mark = [&](nodeid v) {
				if (dont_visit_and_mark_node_if(v))
					return;
				else if (!flc.isReachableFromSource(v)) {
					flc.source_nodes.add_but_dont_assimilate(v);
					visit_node(v);
				}
			};

			auto scan_hyperedge = [&](nodeid from, netid e) {
				if (gather_stats) { bfsstats.scanned_hyperedges++; bfsstats.visited_pins += hg.pinCount(e); }
				assert(is_hyperedge_distance_stale(e));
				current_pin[e] = hg.beginPins(e);
				flc.source_hyperedges.add_but_dont_assimilate(e);
				hyperedge_distance[e] = bfs_dist;
				for (nodeid v : hg.pinsOf(e))
					visit_node_and_mark(v);
			};

			auto visit_hyperedge = [&](nodeid from, netid e) {
				if (gather_stats) { bfsstats.visited_hyperedges++; }
				if (!flc.source_hyperedges.contains(e) && flc.flow_from[e] != from) {
					if (!flc.has_flow(e) || flc.flow_to[e] == from)
						scan_hyperedge(from, e);
					else {
						if (gather_stats) { bfsstats.reroute_attempts++; bfsstats.visited_pins++; }
						visit_node_and_mark(flc.flow_from[e]);
					}
				}
			};

			auto finish_layer = [&]() {
				bfs_dist++;
				queue.finishNextLayer();
			};

			for (nodeid s : flc.source_piercing_nodes)
				visit_node(s);
			finish_layer();

			while (!queue.empty() && !stop_search()) {		//finish scanning current layer in case more hyperedges to target are found.
				while (!queue.currentLayerEmpty()) {
					nodeid u = queue.pop();
					if (gather_stats) { bfsstats.scanned_nodes++; }
					for (netid e : hg.hyperedgesOf(u))
						visit_hyperedge(u, e);
				}
				finish_layer();
			}
			printAndResetBFSStats();
		}

		void growReachableAndBuildLayeredNetworkImpl(const Hypergraph &hg, HLAFlowCutter &flc) {
			auto dont_visit_and_mark_node_if = [](nodeid v) { return false; };
			auto stop_search = []() { return false;};
			buildLayeredNetwork(hg, flc, dont_visit_and_mark_node_if, stop_search);
		}

		bool dinicBFS(const Hypergraph &hg, HLAFlowCutter &flc) {
			reset_distances_if_necessary();
			flc.source_nodes.reset_all_except_assimilated();
			flc.source_hyperedges.reset_all_except_assimilated();
			bool found_target = false;

			auto dont_visit_and_mark_node_if = [&](nodeid v) {
				if (flc.isTarget(v)) {
					node_distance[v] = bfs_dist;
					found_target = true;
					if (gather_stats) { bfsstats.times_seen_target++; }
					return true;
				}
				return false;
			};

			auto stop_search = [&]() { return found_target; };

			buildLayeredNetwork(hg, flc, dont_visit_and_mark_node_if, stop_search);

			return found_target;
		}

		flow_t dinicDFS(const Hypergraph &hg, HLAFlowCutter &flc) {
			flow_t flow = 0;

			auto scan_hyperedge = [&](nodeid u, netid e) {
				if (gather_stats) { dfsstats.hyperedge_scans++; }
				assert(flc.source_hyperedges.contains(e));
				assert(!flc.source_hyperedges.is_assimilated(e));
				assert(!is_hyperedge_distance_stale(e));
				assert(node_distance[u] + 1 == hyperedge_distance[e]);
				for ( ; current_pin[e] != hg.endPins(e); current_pin[e]++) {
					if (gather_stats) { dfsstats.visited_pins++; }
					nodeid v = *current_pin[e];
					assert(flc.isReachableFromSource(v) || flc.isTarget(v));
					if (node_distance[v] == node_distance[u] + 1) {
						assert(flc.isTarget(v) || node_distance[v] == hyperedge_distance[e]);
						return v;
					}
				}
				hyperedge_distance[e] = default_stale_distance;
				return invalid_node;
			};

			for (nodeid s : flc.source_piercing_nodes) {
				while (current_hyperedge[s] != hg.endHyperedges(s)) {
					assert(stack.empty());
					stack.push({s});
					if (gather_stats) { dfsstats.stack_pushes++; }
					while (!stack.empty()) {
						if (gather_stats) { dfsstats.node_scans++; }
						nodeid u = stack.top().tail;
						assert(!is_node_distance_stale(u));
						assert(node_distance[u] - base_distance == stack.size() - 1);
						nodeid v = invalid_node;

						while (current_hyperedge[u] != hg.endHyperedges(u)) {
							netid e = *current_hyperedge[u];
							if (flc.flow_from[e] != u) {
								if (!flc.has_flow(e) || flc.flow_to[e] == u) {
									if (hyperedge_distance[e] == node_distance[u] + 1) {
										v = scan_hyperedge(u, e);
									}
								}
								else {
									if (node_distance[flc.flow_from[e]] == node_distance[u] + 1) {
										if (gather_stats) { dfsstats.reroute_attempts++; }
										v = flc.flow_from[e];
									}
								}
							}

							if (v == invalid_node)
								current_hyperedge[u]++;
							else
								break;
						}

						if (v == invalid_node) {
							if (gather_stats) { dfsstats.stack_pops++; }
							assert(current_hyperedge[u] == hg.endHyperedges(u));
							node_distance[u] = default_stale_distance;
							stack.pop();
						}
						else if (flc.isTarget(v)) {
							//backtrack and push flow.
							if (gather_stats) { dfsstats.flow_diff++; }
							flow++;
							while (!stack.empty()) {
								u = stack.pop().tail;
								assert(current_hyperedge[u] != hg.endHyperedges(u));
								assert(!is_node_distance_stale(u));
								netid e = *current_hyperedge[u];
								Saturate::saturate(flc, u, e, v);
								v = u;
							}
							assert(u == s);
						}
						else {
							if (gather_stats) { dfsstats.stack_pushes++; }
							assert(!is_node_distance_stale(v));
							stack.push( {v} );
						}
					}
				}
			}

			printAndResetDFSStats();
			return flow;
		}

		void swapCutterSpecificDatastructuresAndLocalDatastructures(CutterSpecificDatastructures& csd) {
			std::swap(node_distance, csd.csd_node_distance);
			std::swap(hyperedge_distance, csd.csd_hyperedge_distance);
			std::swap(current_pin, csd.csd_current_pin);
			std::swap(current_hyperedge, csd.csd_current_hyperedge);
		}

		struct BFSStats {
			std::size_t visited_nodes = 0;
			std::size_t scanned_nodes = 0;
			std::size_t visited_hyperedges = 0;
			std::size_t scanned_hyperedges = 0;
			std::size_t visited_pins = 0;
			std::size_t reroute_attempts = 0;
			std::size_t times_seen_target = 0;
		};

		BFSStats bfsstats;

		struct DFSStats {
			std::size_t hyperedge_scans = 0;	//this can be greater than the number of scanned_hyperedges, as they need not be scanned in their entirety from one node. the additional scans incur an additional cache miss compared to BFS
			std::size_t node_scans = 0;
			std::size_t visited_pins = 0;
			std::size_t reroute_attempts = 0;
			std::size_t stack_pushes = 0;
			std::size_t stack_pops = 0;
			flow_t flow_diff = 0;
		};

		DFSStats dfsstats;

		void printBFSStats() {
			if (gather_stats && LoggingInformation::get_output_detail() >= 400) {
				LOG
					<< "BFS "
					<< "vis-n=" << bfsstats.visited_nodes
					<< "scans-n=" << bfsstats.scanned_nodes
					<< "vis-he=" << bfsstats.visited_hyperedges
					<< "scan-he=" << bfsstats.scanned_hyperedges
					<< "try-rr=" << bfsstats.reroute_attempts
					<< "vis-p=" << bfsstats.visited_pins
					<< "see-tar=" << bfsstats.times_seen_target
					<< "tar-depth=" << (bfs_dist - base_distance)
					<< V(bfs_dist) << V(base_distance);
			}
		}

		void printAndResetBFSStats() {
			if (gather_stats) {
				printBFSStats();
				bfsstats = BFSStats();
			}
		}

		void printDFSStats() {
			if (gather_stats && LoggingInformation::get_output_detail() >= 400) {
				LOG
					<< "DFS "
					<< "scan-n=" << dfsstats.node_scans
					<< "scan-he=" << dfsstats.hyperedge_scans
					<< "vis-p=" << dfsstats.visited_pins
					<< "try-rr=" << dfsstats.reroute_attempts
					<< "flow-diff=" << dfsstats.flow_diff;
			}
		}

		void printAndResetDFSStats() {
			if (gather_stats) {
				printDFSStats();
				dfsstats = DFSStats();
			}
		}

		void verify_datastructures_are_released() {
			verify_datastructures_fit_size(0, 0);
		}

		void verify_datastructures_are_acquired() {
			verify_datastructures_fit_size(hg.numNodes(), hg.numHyperedges());
		}

		void verify_datastructures_fit_size(nodeid numNodes, netid numHyperedges) {
			assert(node_distance.size() == numNodes);
			assert(hyperedge_distance.size() == numHyperedges);
			assert(current_hyperedge.size() == numNodes);
			assert(current_pin.size() == numHyperedges);
		}

		void verify_flow_is_exhausted(const Hypergraph& hg, HLAFlowCutter& flc) {
#ifndef NDEBUG
			boost::dynamic_bitset<> src_reachable_nodes(hg.numNodes());
			for (nodeid u = 0; u < hg.numNodes(); u++)
				if (flc.isReachableFromSource(u))
					src_reachable_nodes.set(u);

			boost::dynamic_bitset<> src_scanned_hyperedges(hg.numHyperedges());
			for (netid e = 0; e < hg.numHyperedges(); e++)
				if (flc.source_hyperedges.contains(e))
					src_scanned_hyperedges.set(e);

			flc.flipSearchDirection();
			HLAGrowReachable::growReachable(hg, flc, queue);
			flc.flipSearchDirection();
			HLAGrowReachable::growReachable(hg, flc, queue);

			for (nodeid u = 0; u < hg.numNodes(); u++) {
				assert(!(flc.isReachableFromSource(u) && flc.isTarget(u)));
				assert(src_reachable_nodes[u] == flc.isReachableFromSource(u));
			}
			for (netid e = 0; e < hg.numHyperedges(); e++)
				assert(src_scanned_hyperedges[e] == flc.source_hyperedges.contains(e));
#endif
		}

		void verify_distances_are_stale() {
#ifndef NDEBUG
			for (distance_t dist : node_distance) assert(dist < base_distance);
			for (distance_t dist : hyperedge_distance) assert(dist < base_distance);
#endif
		}


	};



}//namespace hyper
