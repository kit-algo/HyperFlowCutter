#pragma once

#include <tlx/unused.hpp>
#include "hypergraph_search_algorithms.h"
#include "hla_flowcutter.h"
#include "../datastructure/queue.h"
#include "grow_assimilated.h"
#include "saturate.h"


namespace hyper {
	class VertexDisjointHLABFS : public SearchAlgorithm {
	public:
		//static constexpr bool debug = false;
		static constexpr bool debug = true;
		struct LabeledTarget {
			nodeid target;
			label parent_label;
			nodeid parent;
			netid via_hyperedge;
		};

		struct Parent {
			nodeid node;
			netid hyperedge;
		};

		class CutterSpecificDatastructures : public SearchAlgorithm::CutterSpecificDatastructures {
		public:
			std::vector<Parent> parent;
			std::vector<label> plabel;

			explicit CutterSpecificDatastructures(const Hypergraph &hg) : parent(hg.numNodes(), {invalid_node, invalid_hyperedge}),
																		  plabel(hg.numNodes(), invalid_label)
			{

			}

			virtual CutterSpecificDatastructures* clone_impl() const override { return new CutterSpecificDatastructures(*this); }
		};


		bool build_parents_during_grow_reachable;	//ignored here.

		std::vector<LabeledTarget> found_targets;
		LayeredQueue<nodeid> queue;
		std::vector<nodeid> next_layer_label_occurences;

		using Bitvector = boost::dynamic_bitset<>;
		//using Bitvector = Optimization::BoolVec<label>;		//the template instantiation with label is bad if we ever change hyperedge ids to uint64_t but not labels...
		//Optimization::BoolVec shaves off about 8ms from scanning piercing nodes.

		Bitvector label_at_target;    //if the layer containing the first seen target is the last bfs layer, this would be equivalent to next_layer_label_occurences > 0
		Bitvector node_on_next_layer;
		Bitvector hyperedge_reached_target_as_flow_from;

		explicit VertexDisjointHLABFS(const Hypergraph &hg, bool build_datastructures_during_grow_reachable) :
				SearchAlgorithm(hg, build_datastructures_during_grow_reachable),
				build_parents_during_grow_reachable(build_datastructures_during_grow_reachable),
				queue(hg.numNodes()),
				next_layer_label_occurences(std::max(hg.numNodes(), hg.numHyperedges())),
				label_at_target(std::max(hg.numNodes(), hg.numHyperedges())),
				node_on_next_layer(hg.numNodes()),
				hyperedge_reached_target_as_flow_from(hg.numHyperedges()) {}
		//VertexDisjointHLABFS() : VertexDisjointHLABFS(0,0) { }

		std::unique_ptr<SearchAlgorithm::CutterSpecificDatastructures> obtainOwnedCutterSpecificDatastructures() const override {
			return std::make_unique<CutterSpecificDatastructures>(hg);
		}

		void growAssimilated(HLAFlowCutter &flc) override { return HLAAssimilation::growAssimilated(hg, flc, queue); }

		//void growReachableWithoutReset(const Hypergraph& hg, HLAFlowCutter& flc) { return HLAGrowReachable::growReachableWithoutReset(hg, flc, queue); }
		//void growReachable(const Hypergraph& hg, HLAFlowCutter& flc) { return HLAGrowReachable::growReachable(hg, flc, queue); }

		void growReachable(HLAFlowCutter &flc) override {
			flc.source_nodes.reset_all_except_assimilated();
			flc.source_hyperedges.reset_all_except_assimilated();
			growReachableWithoutReset(flc);
		}

		void growReachableWithoutReset(HLAFlowCutter &flc) override {
			growReachableWithoutResetImpl(flc);
			check_datastructures();
		}

		void growReachableWithoutResetImpl(HLAFlowCutter &flc) {
			queue.clear();
			label maxLabel = 0;

			auto & csd = static_cast<CutterSpecificDatastructures&>(*flc.csd);

			auto labelAndPushWithRelabeling = [&](nodeid u, netid e, nodeid v, label lu) {
				assert(!flc.isTarget(v));
				if (!flc.source_nodes.contains(v)) {
					if (!node_on_next_layer[v]) {
						node_on_next_layer.set(v);
						csd.parent[v].node = u;
						csd.parent[v].hyperedge = e;
						csd.plabel[v] = lu;
						next_layer_label_occurences[lu]++;
						queue.push(v);    //not seen yet. push to queue
					} else if (next_layer_label_occurences[csd.plabel[v]] >
							   next_layer_label_occurences[lu] + 1) {    //rebalance the labelings.
						label lv = csd.plabel[v];
						assert(csd.parent[v].node < hg.numNodes());
						assert(csd.parent[v].hyperedge < hg.numHyperedges());
						assert(csd.plabel[v] <= maxLabel);
						assert(next_layer_label_occurences[csd.plabel[v]] > 0);
						csd.parent[v].node = u;
						csd.parent[v].hyperedge = e;
						csd.plabel[v] = lu;
						next_layer_label_occurences[lu]++;    //nllo[lv] does not go below 1! --> packing is tight and no repairs are necessary!
						next_layer_label_occurences[lv]--;
						assert(next_layer_label_occurences[lv] >= 1);
						assert(next_layer_label_occurences[csd.plabel[v]] <= hg.numNodes());    //if due to a bug an underflow occurs, this will catch it
					}
				}
			};

			auto assimilateCurrentLayerAndResetLabelOccurences = [&]() {
				queue.finishNextLayer();
				while (!queue.previousLayerEmpty()) {
					nodeid v = queue.previousLayerPop();
					label lv = csd.plabel[v];
					flc.source_nodes.add_but_dont_assimilate(v);
					if (next_layer_label_occurences[lv] > 0) {
						next_layer_label_occurences[lv] = 0;
					}
				}
			};

			auto bfsTreatHyperedge = [&](nodeid from, netid e, label label_from) {
				if (!flc.source_hyperedges.contains(e) && flc.flow_from[e] != from) {
					if (!flc.has_flow(e) || flc.flow_to[e] == from) {
						flc.source_hyperedges.add_but_dont_assimilate(e);
						for (nodeid v : hg.pinsOf(e)) {
							assert(!flc.isTarget(v));
							labelAndPushWithRelabeling(from, e, v, label_from);
						}
					} else {
						nodeid v = flc.flow_from[e];
						assert(!flc.isTarget(v));
						labelAndPushWithRelabeling(from, e, v, label_from);
					}
				}
			};

			for (nodeid s : flc.source_piercing_nodes) {
				for (netid e : hg.hyperedgesOf(s)) {
					//if (next_layer_label_occurences[maxLabel] > 0) { maxLabel++; }
					maxLabel += static_cast<label>(next_layer_label_occurences[maxLabel] > 0);
					bfsTreatHyperedge(s, e, maxLabel);
				}
			}
			assimilateCurrentLayerAndResetLabelOccurences();

			while (!queue.empty()) {
				while (!queue.currentLayerEmpty()) {
					nodeid u = queue.pop();
					label lu = csd.plabel[u];
					for (netid e : hg.hyperedgesOf(u)) {
						bfsTreatHyperedge(u, e, lu);
					}
				}
				assimilateCurrentLayerAndResetLabelOccurences();
			}
			resetNodesOnNextLayer();
		}

		flow_t extractAugmentingPathsFromPreviousGrowingOfReachableSets(HLAFlowCutter &flc) {
			flow_t flow = 0;
			//don't do it on the first setup iteration
			if (flc.firstGrowReachable || !flc.augmentingPathExists) { return flow; }
			auto& csd = static_cast<CutterSpecificDatastructures&>(*flc.csd);
			std::vector<Parent> path;
			std::vector<label> used_labels;
			assert(!label_at_target.any());
			for (nodeid s : flc.source_piercing_nodes) {
				assert(s < hg.numNodes());
				label ls = csd.plabel[s];
				assert(!flc.isTarget(s));
				assert(ls != invalid_label);
				assert(ls < label_at_target.size());
				if (flc.isReachableFromTarget(s) && !label_at_target[ls]) {
					flow++;
					label_at_target.set(ls);
					used_labels.push_back(ls);
					path.clear();
					nodeid u = s;
					while (!flc.isTarget(u)) {
						assert(csd.plabel[u] == ls);
						assert(u < hg.numNodes());
						assert(u < csd.parent.size());
						nodeid v = csd.parent[u].node;
						path.push_back(csd.parent[u]);
						u = v;
					}
					assert(u == path.back().node);
					nodeid target = u;
					tlx::unused(target);    //only for assertion
					u = s;
					for (auto p : path) {
						Saturate::saturate(flc, u, p.hyperedge, p.node);
						u = p.node;
					}
					assert(u == target);
				}
			}
			for (label l : used_labels) { label_at_target.reset(l); }
			check_datastructures();
			return flow;
		}

		flow_t exhaustFlow(HLAFlowCutter &flc) override {
			auto & csd = static_cast<CutterSpecificDatastructures&>(*flc.csd);
			flow_t flow = 0;
			flow += extractAugmentingPathsFromPreviousGrowingOfReachableSets(flc);
			while (growFlows(flc)) {
				auto flow_diff = found_targets.size();
				flow += flow_diff;
				for (const auto &x : found_targets) {
					assert(label_at_target[x.parent_label]);
					assert(flc.isTarget(x.target));
					label_at_target.reset(x.parent_label);
					next_layer_label_occurences[x.parent_label] = 0;
					hyperedge_reached_target_as_flow_from.reset(x.via_hyperedge);
					nodeid v = x.target;
					nodeid u = x.parent;
					netid e = x.via_hyperedge;
					Saturate::saturate(flc, u, e,
									   v);    //saturate the hyperedge entering target separately since it is not in the parent pointer
					v = u;
					while (!flc.isSource(v)) {
						u = csd.parent[v].node;
						e = csd.parent[v].hyperedge;
						assert(csd.plabel[u] == x.parent_label || flc.isSource(u));
						assert(u != invalid_node);
						assert(e != invalid_hyperedge);
						Saturate::saturate(flc, u, e, v);
						v = u;
					}
				}
				found_targets.clear();
				assert(!label_at_target.any());
				check_datastructures();
			}
			return flow;
		}


		bool growFlows(HLAFlowCutter &flc) {
			auto & csd = static_cast<CutterSpecificDatastructures&>(*flc.csd);
			found_targets.clear();
			flc.source_nodes.reset_all_except_assimilated();
			flc.source_hyperedges.reset_all_except_assimilated();
			queue.clear();
			label numberOfActiveLabels = 0;
			label maxLabel = 0;

			auto labelAndPushWithRelabeling = [&](nodeid u, netid e, nodeid v, label lu) {
				assert(!label_at_target[lu]);
				assert(!flc.isTarget(v));
				if (!flc.source_nodes.contains(v)) {
					if (!node_on_next_layer[v]) {
						node_on_next_layer.set(v);
						csd.parent[v].node = u;
						csd.parent[v].hyperedge = e;
						csd.plabel[v] = lu;
						next_layer_label_occurences[lu]++;
						queue.push(v);    //not seen yet. push to queue
					} else if (label_at_target[csd.plabel[v]] || next_layer_label_occurences[csd.plabel[v]] >
																 next_layer_label_occurences[lu] +
																 1) {    //old label reaches target -> relabel. or rebalance the labelings.
						label lv = csd.plabel[v];
						assert(csd.parent[v].node < hg.numNodes());
						assert(csd.parent[v].hyperedge < hg.numHyperedges());
						assert(csd.plabel[v] <= maxLabel);
						assert(next_layer_label_occurences[csd.plabel[v]] > 0);
						csd.parent[v].node = u;
						csd.parent[v].hyperedge = e;
						csd.plabel[v] = lu;
						next_layer_label_occurences[lu]++;    //nllo[lv] does not go below 1! --> packing is tight and no repairs are necessary!
						next_layer_label_occurences[lv]--;
						assert(next_layer_label_occurences[lv] >= 1);
						assert(next_layer_label_occurences[csd.plabel[v]] <=
							   hg.numNodes());    //if due to a bug an underflow occurs, this will catch it
					}
				}
			};

			auto isTargetLabelUnusedAndHyperedgeUsable = [&](nodeid u, netid e, nodeid v, label lu) {
				assert(flc.isTarget(v));
				if (!label_at_target[lu] && !(flc.flow_from[e] == v && hyperedge_reached_target_as_flow_from[e])) {
					label_at_target.set(lu);
					next_layer_label_occurences[lu]++;
					found_targets.push_back({v, lu, u, e});
					if (flc.flow_from[e] == v) { hyperedge_reached_target_as_flow_from.set(e); }
					return true;
				}
				return false;
			};

			auto assimilateCurrentLayerResetLabelOccurencesAndComputeNumberOfActiveLabels = [&]() {
				numberOfActiveLabels = 0;
				queue.finishNextLayer();
				while (!queue.previousLayerEmpty()) {
					nodeid v = queue.previousLayerPop();
					label lv = csd.plabel[v];
					flc.source_nodes.add_but_dont_assimilate(v);
					if (next_layer_label_occurences[lv] > 0) {
						next_layer_label_occurences[lv] = 0;
						if (!label_at_target[lv]) { numberOfActiveLabels++; }
					}
				}
			};

			auto bfsTreatHyperedge = [&](nodeid from, netid e, label label_from) {
				if (!flc.source_hyperedges.contains(e) && flc.flow_from[e] != from) {
					if (!flc.has_flow(e) || flc.flow_to[e] == from) {
						flc.source_hyperedges.add_but_dont_assimilate(e);
						for (nodeid v : hg.pinsOf(e)) {
							if (flc.isTarget(v)) {
								if (isTargetLabelUnusedAndHyperedgeUsable(from, e, v, label_from)) { return false; }
							}
							else { labelAndPushWithRelabeling(from, e, v, label_from); }
						}
					} else {
						nodeid v = flc.flow_from[e];
						if (flc.isTarget(v)) {
							if (isTargetLabelUnusedAndHyperedgeUsable(from, e, v, label_from)) { return false; }
						}
						else { labelAndPushWithRelabeling(from, e, v, label_from); }
					}
				}
				return true;
			};

			{
				for (nodeid s : flc.source_piercing_nodes) {
					for (netid e : hg.hyperedgesOf(s)) {
						maxLabel += static_cast<label>(next_layer_label_occurences[maxLabel] > 0);
						bfsTreatHyperedge(s, e, maxLabel);
					}
				}
				assimilateCurrentLayerResetLabelOccurencesAndComputeNumberOfActiveLabels();
			}

			{
				while (!queue.empty() && numberOfActiveLabels > 0) {
					while (!queue.currentLayerEmpty()) {
						nodeid u = queue.pop();
						label lu = csd.plabel[u];
						if (label_at_target[lu]) { continue; }
						for (netid e : hg.hyperedgesOf(u)) {
							if (!bfsTreatHyperedge(u, e, lu))
								break;
						}
					}
					assimilateCurrentLayerResetLabelOccurencesAndComputeNumberOfActiveLabels();
				}
			}
			resetNodesOnNextLayer();
			assert(label_at_target.count() == found_targets.size());
			return !found_targets.empty();
		}

		inline void resetNodesOnNextLayer() {
			if (queue.numberOfPushedElements() > node_on_next_layer.size() / Bitvector::bits_per_block) {
				node_on_next_layer.reset();
			} else {
				queue.forAllEverContainedElements([&](nodeid u) { node_on_next_layer.reset(u); });
			}
		}

		void check_datastructures() {
#ifndef NDEBUG
			for (auto x : next_layer_label_occurences) {
				assert(x == 0);
			}
			assert(!node_on_next_layer.any());
			assert(!label_at_target.any());
			assert(!hyperedge_reached_target_as_flow_from.any());
			assert(found_targets.empty());
#endif
		}

	};

}//namespace