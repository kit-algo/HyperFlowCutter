#pragma once

#include <unordered_map>
#include <cmath>
#include <tlx/unused.hpp>
#include "../definitions.h"
#include "../datastructure/mixed_assimilation_timestamp_set.h"
#include "boost/dynamic_bitset.hpp"
#include "../util/timer.h"
#include "../datastructure/hypergraph.h"
#include "../datastructure/sparse_reset_vector.h"
#include "../datastructure/partition.h"
#include "hypergraph_search_algorithms.h"


namespace hyper {

	class HLAFlowCutter {
	public:

		using timestamp_t = uint8_t;
		using search_space= assimilation_timestamp_set<nodeid, timestamp_t>;
		using hyperedge_search_space = hyperedge_assimilation_timestamp_set<netid, timestamp_t>;
		using Bitvector = boost::dynamic_bitset<>;

		//using search_space = BitVectorAssimilationSet<nodeid>;
		//using hyperedge_search_space = BitVectorAssimilationSet<netid>;

		int step = 0;

		static constexpr int source_side = 0;
		static constexpr int target_side = 1;

		static constexpr bool debug = false;

		bool augmentingPathExists = true;
		int piercingSide = -1;

		bool firstGrowReachable = false;

		nodeid total_node_weight = 0, half_total_node_weight_ceiled = 0;
		nodeid contracted_source_assimilated_nodes = 0, contracted_target_assimilated_nodes = 0;

		nodeid desired_smaller_blocksize = 0;
		bool hit_desired_smaller_blocksize_exactly = false;

		int source_side_in_current_direction = 0;
		flow_t cutsize = 0;

		std::vector<nodeid> flow_from, flow_to;
		//Bitvector has_flow;	//TODO can be eliminated by using flow_from != invalid_node



		search_space source_nodes, target_nodes;						//could be turned into one since reachable sets must be disjoint, unless the node got pierced. need separate timestamps for source_assimilated, target_assimilated
		hyperedge_search_space source_hyperedges, target_hyperedges;	//could be turned into one. odd timestamps for source, even for target
		std::vector<nodeid> source_piercing_nodes, target_piercing_nodes;
		Bitvector was_added_to_source_cut, was_added_to_target_cut;
		Bitvector has_source_assimilated_pins, has_target_assimilated_pins;


		std::vector<nodeid> active_pincount;
		//std::vector<netid> number_of_hyperedges_with_source_pins, number_of_hyperedges_with_target_pins;

		std::vector<netid> source_cut_front, target_cut_front;
		std::vector<netid> source_mixed_cut, target_mixed_cut;
		std::vector<nodeid> source_border_vertices, target_border_vertices;
		Bitvector was_added_to_source_border, was_added_to_target_border;

		nodeid nIsolatedNodes = 0;
		std::vector<netid> number_of_mixed_incident_hyperedges;
		std::vector<nodeid> free_nodes;
		std::vector<nodeid> free_nodes_candidates;
		OverlapWithPartitionBlocks source_side_overlap_with_ensemble_partition, target_side_overlap_with_ensemble_partition;

		hyper::duration running_time = duration(0.0);
		SearchAlgorithm& search_algorithm;
		std::unique_ptr<SearchAlgorithm::CutterSpecificDatastructures> csd;
		const Hypergraph& hg;

		//HLAFlowCutter() = default;

		explicit HLAFlowCutter(const Hypergraph& hg, SearchAlgorithm& search_algorithm, OverlapWithPartitionBlocks& ensemble_overlap, bool hit_desired_smaller_blocksize_exactly = false) :
				augmentingPathExists(true),
				total_node_weight(hg.totalNodeWeight()),
				half_total_node_weight_ceiled(hg.totalNodeWeight() / 2 + (hg.totalNodeWeight() % 2 == 0 ? 0 : 1)  ),
				desired_smaller_blocksize(hg.totalNodeWeight() / 2),
				hit_desired_smaller_blocksize_exactly(hit_desired_smaller_blocksize_exactly),
				source_side_in_current_direction(source_side),
				cutsize(0),
				flow_from(hg.numHyperedges(), invalid_node), flow_to(hg.numHyperedges(), invalid_node),
				source_nodes(hg.numNodes()), target_nodes(hg.numNodes()),
				source_hyperedges(hg.numHyperedges()), target_hyperedges(hg.numHyperedges()),
				was_added_to_source_cut(hg.numHyperedges()), was_added_to_target_cut(hg.numHyperedges()),
				has_source_assimilated_pins(hg.numHyperedges()), has_target_assimilated_pins(hg.numHyperedges()),
				active_pincount(hg.numHyperedges()),
				was_added_to_source_border(hg.numNodes()), was_added_to_target_border(hg.numNodes()),
				number_of_mixed_incident_hyperedges(hg.numNodes(), 0),
				source_side_overlap_with_ensemble_partition(ensemble_overlap), target_side_overlap_with_ensemble_partition(ensemble_overlap),
				search_algorithm(search_algorithm),
				csd(search_algorithm.obtainOwnedCutterSpecificDatastructures()),
				hg(hg)
		{
			for (netid e = 0; e < hg.numHyperedges(); e++) { active_pincount[e] = hg.pinCount(e); }
		}

		HLAFlowCutter(const HLAFlowCutter& other) : step(other.step), augmentingPathExists(other.augmentingPathExists), piercingSide(other.piercingSide), firstGrowReachable(other.firstGrowReachable), total_node_weight(other.total_node_weight),
													half_total_node_weight_ceiled(other.half_total_node_weight_ceiled), contracted_source_assimilated_nodes(other.contracted_source_assimilated_nodes), contracted_target_assimilated_nodes(other.contracted_target_assimilated_nodes),
													desired_smaller_blocksize(other.desired_smaller_blocksize), hit_desired_smaller_blocksize_exactly(other.hit_desired_smaller_blocksize_exactly), source_side_in_current_direction(other.source_side_in_current_direction),
													cutsize(other.cutsize), flow_from(other.flow_from), flow_to(other.flow_to), source_nodes(other.source_nodes), target_nodes(other.target_nodes), source_hyperedges(other.source_hyperedges), target_hyperedges(other.target_hyperedges),
													source_piercing_nodes(other.source_piercing_nodes), target_piercing_nodes(other.target_piercing_nodes), was_added_to_source_cut(other.was_added_to_source_cut), was_added_to_target_cut(other.was_added_to_target_cut),
													has_source_assimilated_pins(other.has_source_assimilated_pins), has_target_assimilated_pins(other.has_target_assimilated_pins), active_pincount(other.active_pincount), source_cut_front(other.source_cut_front), target_cut_front(other.target_cut_front),
													source_mixed_cut(other.source_mixed_cut), target_mixed_cut(other.target_mixed_cut), source_border_vertices(other.source_border_vertices), target_border_vertices(other.target_border_vertices), was_added_to_source_border(other.was_added_to_source_border),
													was_added_to_target_border(other.was_added_to_target_border), nIsolatedNodes(other.nIsolatedNodes), number_of_mixed_incident_hyperedges(other.number_of_mixed_incident_hyperedges), free_nodes(other.free_nodes),
													free_nodes_candidates(other.free_nodes_candidates), source_side_overlap_with_ensemble_partition(other.source_side_overlap_with_ensemble_partition), target_side_overlap_with_ensemble_partition(other.target_side_overlap_with_ensemble_partition),
													running_time(other.running_time),
													search_algorithm(other.search_algorithm),
													csd(other.csd->clone()),		//this single line is the only part where something happens. the rest of this constructor is nasty error-prone boilerplate
													hg(other.hg)
		{	}


		inline void addToCut(netid e) {
			assert(!wasAddedToCut(e));
			for (nodeid pin : hg.pinsOf(e)) {
				if (!isSource(pin) && !isTarget(pin) && !was_added_to_source_border[pin]) {
					source_border_vertices.push_back(pin);
					was_added_to_source_border.set(pin);
				}
			}
			was_added_to_source_cut.set(e);
			source_cut_front.push_back(e);
		}

		inline bool has_flow(const netid e) const { return flow_from[e] != invalid_node; }
		
		inline bool isFreeNode(nodeid u) const { return number_of_mixed_incident_hyperedges[u] == hg.degree(u) && !isSource(u) && !isTarget(u); }

		inline bool isCutHyperedgeReachableFromTarget(netid e) const {
			assert(has_flow(e));
			assert(flow_to[e] != invalid_node);
#ifndef NDEBUG
			bool expected = isReachableFromTarget(flow_to[e]);
			bool has_reachable_pin = false;
			for (nodeid pin : hg.pinsOf(e)) { has_reachable_pin |= isReachableFromTarget(pin);}
			assert(expected == has_reachable_pin);
#endif
            return isReachableFromTarget(flow_to[e]);	//flow_to[e] for the current source-side is flow_from[e] for when the current target-side is source-side --> (any pin of e reachable ==> flow_to[e] is reachable).
        }

		inline bool finishSourceSideDirectly(nodeid blocksize) const {
			return numberOfSourceReachableNodes() <= half_total_node_weight_ceiled && numberOfSourceReachableNodes() + nIsolatedNodes >= blocksize
					&& (!hit_desired_smaller_blocksize_exactly || numberOfSourceReachableNodes() <= blocksize);
		}
		inline bool finishTargetSideDirectly(nodeid blocksize) const {
			return numberOfTargetReachableNodes() <= half_total_node_weight_ceiled && numberOfTargetReachableNodes() + nIsolatedNodes >= blocksize
					&& (!hit_desired_smaller_blocksize_exactly || numberOfTargetReachableNodes() <= blocksize);
		}
		inline bool finishSourceSideWithUnreachable(nodeid blocksize) const {
			return numberOfSourceReachableNodes() + numberOfUnreachableNodes() >= blocksize && numberOfSourceReachableNodes() + numberOfUnreachableNodes() - nIsolatedNodes <= half_total_node_weight_ceiled
				   && (!hit_desired_smaller_blocksize_exactly || numberOfSourceReachableNodes() + numberOfUnreachableNodes() - nIsolatedNodes <= blocksize);
		}
		inline bool finishTargetSideWithUnreachable(nodeid blocksize) const {
			return numberOfTargetReachableNodes() + numberOfUnreachableNodes() >= blocksize && numberOfTargetReachableNodes() + numberOfUnreachableNodes() - nIsolatedNodes <= half_total_node_weight_ceiled
				   && (!hit_desired_smaller_blocksize_exactly || numberOfTargetReachableNodes() + numberOfUnreachableNodes() - nIsolatedNodes <= blocksize);
		}
		inline bool isBlocksizeAchieved(nodeid blocksize) const {
			assert(totalNodeWeight() >= numberOfSourceReachableNodes() + numberOfTargetReachableNodes());
			return finishSourceSideDirectly(blocksize) || finishTargetSideDirectly(blocksize) || finishSourceSideWithUnreachable(blocksize) || finishTargetSideWithUnreachable(blocksize);
		}

		inline bool isFullyAssimilated() const { return numberOfSourceAssimilatedNodes() + numberOfTargetAssimilatedNodes() + nIsolatedNodes == totalNodeWeight(); }

		inline bool isDesiredlyBalanced() const { return isBlocksizeAchieved(desired_smaller_blocksize); }
		inline bool isEpsilonBalanced(double epsilon) const {
			return isBlocksizeAchieved( Metrics::smallerBlockSize(totalNodeWeight(), epsilon) );
		}

		inline double currentImbalance() const {
			if (isBlocksizeAchieved(totalNodeWeight() / 2)) { return 0.0; }
			nodeid sma = 0;
			//maximize the smaller block!
			if (numberOfSourceReachableNodes() <= totalNodeWeight()/2) { sma = std::max(sma, std::min(totalNodeWeight()/2, numberOfSourceReachableNodes() + nIsolatedNodes)); }
			if (numberOfTargetReachableNodes() <= totalNodeWeight()/2) { sma = std::max(sma, std::min(totalNodeWeight()/2, numberOfTargetReachableNodes() + nIsolatedNodes)); }
			if (numberOfSourceReachableNodes() + numberOfUnreachableNodes() - nIsolatedNodes <= totalNodeWeight()/2) {
				sma = std::max(sma, std::min(totalNodeWeight()/2, numberOfSourceReachableNodes() + numberOfUnreachableNodes()));
			}
			if (numberOfTargetReachableNodes() + numberOfUnreachableNodes() - nIsolatedNodes <= totalNodeWeight()/2) {
				sma = std::max(sma, std::min(totalNodeWeight()/2, numberOfTargetReachableNodes() + numberOfUnreachableNodes()));
			}

			auto smallerside = static_cast<double>(sma);
			return 1.0 - 2.0 * smallerside / static_cast<double>(totalNodeWeight());
		}

		int sideToPierce() const {
			nodeid num_nodes_smr = std::min(numberOfSourceReachableNodes(), numberOfTargetReachableNodes());
			if (hit_desired_smaller_blocksize_exactly && num_nodes_smr > desired_smaller_blocksize) {
				return largerAssimilatedSide();	//at this stage, assimilating either side would break that side's opportunity to be the smaller side. so if one side is already grown beyond desired_smaller_blocksize, we must grow that side again.
			}
			else {
				return smallerReachableSide();
			}
		}

		int largerAssimilatedSide() const { return numberOfSourceAssimilatedNodes() >= numberOfTargetAssimilatedNodes() ? source_side_in_current_direction : 1 - source_side_in_current_direction; }
		int smallerReachableSide() const { return numberOfSourceReachableNodes() <= numberOfTargetReachableNodes() ? source_side_in_current_direction : 1-source_side_in_current_direction; }
		int cutSide() const {
			//source_side fully assimilated
			if (numberOfSourceAssimilatedNodes() == numberOfSourceReachableNodes()
				&& (
						//and either target side not
						numberOfTargetAssimilatedNodes() != numberOfTargetReachableNodes()
						//or if both are fully assimilated, choose the smaller
						|| numberOfSourceAssimilatedNodes() <= numberOfTargetAssimilatedNodes())
					)
				return source_side_in_current_direction;
			else	//target_side is fully assimilated (as at least one of the two sides is always fully assimilated). AND if both are fully assimilated, target side is smaller
				return 1-source_side_in_current_direction;
		}
		void setSearchDirection(int source_side_in_new_direction) {
			if (source_side_in_current_direction != source_side_in_new_direction) {
				std::swap(contracted_source_assimilated_nodes, contracted_target_assimilated_nodes);
				flow_from.swap(flow_to);
				source_nodes.swap(target_nodes);
				source_hyperedges.swap(target_hyperedges);
				source_cut_front.swap(target_cut_front);
				source_mixed_cut.swap(target_mixed_cut);
				source_border_vertices.swap(target_border_vertices);
				was_added_to_source_border.swap(was_added_to_target_border);
				//number_of_hyperedges_with_source_pins.swap(number_of_hyperedges_with_target_pins);
				source_piercing_nodes.swap(target_piercing_nodes);
				was_added_to_source_cut.swap(was_added_to_target_cut);
				has_source_assimilated_pins.swap(has_target_assimilated_pins);
				source_side_overlap_with_ensemble_partition.swap(target_side_overlap_with_ensemble_partition);
				source_side_in_current_direction = source_side_in_new_direction;
			}
		}
		void flipSearchDirection() { setSearchDirection(1-source_side_in_current_direction); }

		inline bool hasSourceAssimilatedPin(netid e) const { return has_source_assimilated_pins[e]; }
		inline bool hasTargetAssimilatedPin(netid e) const { return has_target_assimilated_pins[e]; }
		inline bool hasAssimilatedPinsOnBothSides(netid e) const { return hasSourceAssimilatedPin(e) && hasTargetAssimilatedPin(e); }

		inline void node_assimilate(nodeid u) {
			assert(!isTarget(u)); assert(!isSource(u));
			source_nodes.assimilate(u);
			if (source_side_overlap_with_ensemble_partition.use_ensemble_piercing) { source_side_overlap_with_ensemble_partition.settleNode(u); }
			for (netid e : hg.hyperedgesOf(u)) {
				active_pincount[e]--;
				if (!hasSourceAssimilatedPin(e)) {
					has_source_assimilated_pins.set(e);
					//for (nodeid pin : hg.pinsOf(e)) { number_of_hyperedges_with_source_pins[pin]++; }
					if (hasTargetAssimilatedPin(e)) {
						for (nodeid pin : hg.pinsOf(e)) {
							number_of_mixed_incident_hyperedges[pin]++;
							if (!isSource(pin) && !isTarget(pin) && number_of_mixed_incident_hyperedges[pin] == hg.degree(pin)) {
								free_nodes_candidates.push_back(pin);
							}
						}
					}
				}
			}
		}

		inline bool shouldBeAddedToCutfront(netid e) const { return !source_hyperedges.contains(e) && !wasAddedToCut(e) && has_flow(e); }

		void deleteNonCutHyperedges() {
			for (netid i = 0; i < source_cut_front.size(); ) {
				assert(source_hyperedges.is_assimilated(source_cut_front[i]) == (active_pincount[source_cut_front[i]] == 0 && !hasAssimilatedPinsOnBothSides(source_cut_front[i])));
				if (hasAssimilatedPinsOnBothSides(source_cut_front[i])) {
					source_mixed_cut.push_back(source_cut_front[i]);
					source_cut_front[i] = source_cut_front.back(); source_cut_front.pop_back();
				}
				else { i++; }
			}
			source_cut_front.erase(
					std::remove_if(source_cut_front.begin(), source_cut_front.end(),
								   [&](const netid e) {
										return source_hyperedges.is_assimilated(e) || active_pincount[e] == 0;
								   }),
					source_cut_front.end());
		}

		void computeFreeNodes() {
			for (nodeid u : free_nodes_candidates) {
				if (!isSource(u) && !isTarget(u)) { free_nodes.push_back(u); nIsolatedNodes++; }
			}
			free_nodes_candidates.clear();
		}


		void growReachable() {
			assert(numberOfSourceAssimilatedNodes() == numberOfSourceReachableNodes());
			if (firstGrowReachable || augmentingPathExists) {
				cutsize += search_algorithm.exhaustFlow(*this);
				flipSearchDirection();
				search_algorithm.growReachable(*this);
			}
			else {
				//no point in growing T_r because a vertex not reachable from T was added to S. --> only grow S_r.
				search_algorithm.growReachableWithoutReset(*this);
			}
			firstGrowReachable = false;
			check_reachable_sides_are_disjoint();
		}

		void growAssimilated() {
			int side = sideToPierce();//smallerReachableSide();
			setSearchDirection(side);
			search_algorithm.growAssimilated(*this);
			deleteNonCutHyperedges();
			computeFreeNodes();
			assert(static_cast<std::size_t>(cutsize) == source_cut_front.size() + source_mixed_cut.size());
			assert(numberOfSourceAssimilatedNodes() == numberOfSourceReachableNodes());
			check_every_cut_hyperedge_has_source_pin();
			check_free_nodes_are_not_reachable();
		}

		void prepareToPierceSingleNodes() {
			source_border_vertices.erase(
					std::remove_if(source_border_vertices.begin(), source_border_vertices.end(), [&](const nodeid pin) { return isSource(pin) || isTarget(pin) || isFreeNode(pin); }),
					source_border_vertices.end() );
			check_border_vertices_are_not_assimilated();
		}

		void pierceHyperedge(netid piercing_hyperedge) {
			hyperedge_assimilate(piercing_hyperedge);
			for (nodeid u : hg.pinsOf(piercing_hyperedge)) {
				if (!isSource(u)) {
					node_assimilate(u);
					source_piercing_nodes.push_back(u);
				}
			}
		}


		template<typename PierceScorer>
		bool findSomethingToPierce(PierceScorer& scorer) {
			source_piercing_nodes.clear();
			netid piercing_hyperedge = invalid_hyperedge;
			nodeid p;

			auto [eligible_piercing_hyperedge_exists, foundPiercingHyperedgeAvoidingAugmentingPaths] = scorer.findAvoidAugmentingPathPiercingHyperedge(*this, hg, piercing_hyperedge);
			tlx::unused(eligible_piercing_hyperedge_exists);
			if (foundPiercingHyperedgeAvoidingAugmentingPaths) {
				//augmentingPathExists = false; welp. this is a bad idea if we allow to ignore aap information...
				augmentingPathExists = isCutHyperedgeReachableFromTarget(piercing_hyperedge);
				pierceHyperedge(piercing_hyperedge); return true;
			}
			else {
				prepareToPierceSingleNodes();
				p = scorer.findAvoidAugmentingPathPiercingNode(*this, hg);
				if (p != invalid_node) {
					source_piercing_nodes.push_back(p);
					node_assimilate(p);
					augmentingPathExists = isReachableFromTarget(p);
					return true;
				}
				else {
					bool eligible_piercing_nodes_exist = !source_border_vertices.empty();
					if (scorer.findAnyPiercingHyperedge(*this, hg, piercing_hyperedge)) {
						pierceHyperedge(piercing_hyperedge);
						augmentingPathExists = isCutHyperedgeReachableFromTarget(piercing_hyperedge);
						return true;
					}
					else if (eligible_piercing_nodes_exist) {
						p = scorer.findAnyPiercingNode(*this, hg);
						if (p != invalid_node) {
							node_assimilate(p);
							source_piercing_nodes.push_back(p);
							augmentingPathExists = isReachableFromTarget(p);
							return true;
						}
					}
				}
			}

			{
				//Potential other option would be to pierce the other side if it has fewer than half of the vertices.
				LOG << "HLAFC: No piercing node from the border possible.";
				//If not desiredly balanced, choose a random non-terminal vertex to assimilate to s. Preferably unreachable from t...but that's not really possible
				if (!isDesiredlyBalanced() && !isFullyAssimilated()) {
					std::vector<nodeid> reachable_candidates, unreachable_candidates;
					for (nodeid u = 0; u < hg.numNodes(); u++)
						if (!isSource(u) && !isTarget(u) && !isFreeNode(u)) {
							if (isReachableFromTarget(u)) reachable_candidates.push_back(u);
							else unreachable_candidates.push_back(u);
						}
					if (unreachable_candidates.empty()) {
						if (reachable_candidates.empty()) {
							LOG << "HLAFC: Even finding a random vertex to pierce has failed!";
							return false;
						}
						std::uniform_int_distribution<nodeid> dist(0, static_cast<nodeid>(reachable_candidates.size() - 1));
						p = reachable_candidates[dist(Random::getRNG())];
					}
					else {
						std::uniform_int_distribution<nodeid> dist(0, static_cast<nodeid>(unreachable_candidates.size()-1));
						p = unreachable_candidates[dist(Random::getRNG())];
					}
					augmentingPathExists = isReachableFromTarget(p);
					node_assimilate(p);
					source_piercing_nodes.push_back(p);
					return true;
				}
			}
			return false;
		}

		void initializeSourceTarget(nodeid s, nodeid t) {
			source_side_in_current_direction = target_side;
			node_assimilate(t);
			flipSearchDirection();
			node_assimilate(s);
			source_piercing_nodes={s}; target_piercing_nodes={t};
		}

		void assimilateTerminalsButOnlyAddToPiercingNodesIfReachableFromBorder(const std::vector<nodeid> &terminals) {
			for (nodeid tx : terminals)
				if (!isSource(tx))
					node_assimilate(tx);
			for (netid e = 0; e < hg.numHyperedges(); e++)
				if (activePinCount(e) == 0 && hasSourceAssimilatedPin(e) && !hasTargetAssimilatedPin(e))
					hyperedge_assimilate(e);
			for (nodeid tx : terminals)
				if (std::any_of(hg.hyperedgesOf(tx).begin(), hg.hyperedgesOf(tx).end(), [&](const netid& e) { return !source_hyperedges.is_assimilated(e); }))
					source_piercing_nodes.push_back(tx);
		}

		void initializeSourceTarget(const std::vector<nodeid>& s, const std::vector<nodeid>& t) {
			source_piercing_nodes.clear();
			target_piercing_nodes.clear();
			source_side_in_current_direction = target_side;
			assimilateTerminalsButOnlyAddToPiercingNodesIfReachableFromBorder(t);
			flipSearchDirection();
			assimilateTerminalsButOnlyAddToPiercingNodesIfReachableFromBorder(s);
		}

		void initialCut(const STPair& stp) {
			if (stp.s.empty() || stp.t.empty()) throw std::runtime_error("Terminal set is empty");
			if ((stp.s_is_coarse && stp.s.size() != 1) || (stp.t_is_coarse && stp.t.size() != 1)) throw std::runtime_error("Coarse terminal has more than one node");
			if (stp.s_is_coarse)
				contracted_source_assimilated_nodes = hg.nodeWeight(stp.s[0]) - 1;
			if (stp.t_is_coarse)
				contracted_target_assimilated_nodes = hg.nodeWeight(stp.t[0]) - 1;

			initializeSourceTarget(stp.s, stp.t);
			cutsize = 0;
			check_fc_invariants();
			growReachable();
			check_fc_invariants();
			growAssimilated();
			check_fc_invariants();
		}

		template<typename PierceScorer>
		bool advance(PierceScorer& scorer) {
			if (isDesiredlyBalanced() || isFullyAssimilated())
				return false;
			check_fc_invariants();
			piercingSide = cutSide();
			setSearchDirection(piercingSide);
			check_fc_invariants();
			if (!findSomethingToPierce(scorer))
				return false;
			assert(augmentingPathExists == std::any_of(source_piercing_nodes.begin(), source_piercing_nodes.end(), [&](const nodeid& src_p) { return isReachableFromTarget(src_p); }));
			++step;
			check_fc_invariants();
			growReachable();
			check_fc_invariants();
			growAssimilated();
			check_fc_invariants();
			return true;
		}

		inline nodeid totalNodeWeight() const { return total_node_weight; }
		inline nodeid numberOfContractedSourceAssimilatedNodes() const { return contracted_source_assimilated_nodes; }
		inline nodeid numberOfContractedTargetAssimilatedNodes() const { return contracted_target_assimilated_nodes; }
		inline nodeid numberOfSourceReachableNodes() const { return source_nodes.num_elements() + numberOfContractedSourceAssimilatedNodes(); }
		inline nodeid numberOfSourceAssimilatedNodes() const { return source_nodes.num_assimilated_elements() + numberOfContractedSourceAssimilatedNodes(); }
		inline nodeid numberOfTargetReachableNodes() const { return target_nodes.num_elements() + numberOfContractedTargetAssimilatedNodes(); }
		inline nodeid numberOfTargetAssimilatedNodes() const { return target_nodes.num_assimilated_elements() + numberOfContractedTargetAssimilatedNodes(); }
		inline nodeid numberOfUnreachableNodes() const {
			assert(totalNodeWeight() >= numberOfSourceReachableNodes() + numberOfTargetReachableNodes());
			return totalNodeWeight() - numberOfSourceReachableNodes() - numberOfTargetReachableNodes();
		}
		inline bool isSource(nodeid u) const { return source_nodes.is_assimilated(u); }
		inline bool isTarget(nodeid u) const { return target_nodes.is_assimilated(u); }
		inline bool isReachableFromSource(nodeid u) const { return source_nodes.contains(u); }
		inline bool isReachableFromTarget(nodeid u) const { return target_nodes.contains(u); }
		inline bool wasAddedToCut(netid e) const { return was_added_to_source_cut[e]; }
		inline void hyperedge_assimilate(netid e) { source_hyperedges.assimilate(e); }
		inline const search_space& nodeSearchSpace(int side) const { if (side == source_side_in_current_direction) return source_nodes; else return target_nodes; }
		inline const hyperedge_search_space& hyperedgeSearchSpace(int side) const { if (side == source_side_in_current_direction) return source_hyperedges; else return target_hyperedges; }
		inline nodeid activePinCount(netid e) const { return active_pincount[e]; }

		inline nodeid numNodes() const { return source_nodes.capacity(); }
		inline netid numHyperedges() const { return source_hyperedges.capacity(); }

		std::string state() {
			nodeid sr = numberOfSourceReachableNodes(), sa = numberOfSourceAssimilatedNodes(), tr = numberOfTargetReachableNodes(), ta = numberOfTargetAssimilatedNodes();
			nodeid scr = source_nodes.num_elements(), sca = source_nodes.num_assimilated_elements(), tcr = target_nodes.num_elements(), tca = target_nodes.num_assimilated_elements();
			auto sPierce = static_cast<nodeid>(source_piercing_nodes.size()), tPierce = static_cast<nodeid>(target_piercing_nodes.size());
			if (source_side_in_current_direction == 1) {
				std::swap(sr, tr); std::swap(sa, ta); std::swap(sPierce, tPierce);
				std::swap(scr, tcr); std::swap(sca, tca);
			}
			return
					"n=" + std::to_string(totalNodeWeight())
					+ " cut=" + std::to_string(cutsize)
					+ " sC=" + std::to_string(sca) + "|" + std::to_string(scr)
					+ " tC=" + std::to_string(tca) + "|" + std::to_string(tcr)
					+ " s: " + std::to_string(sa) + "|" + std::to_string(sr)
					+ " t: " + std::to_string(ta) + "|" + std::to_string(tr)
					+ " nIso=" + std::to_string(nIsolatedNodes)
					+ " sPierce=" + std::to_string(sPierce) + " | tPierce=" + std::to_string(tPierce);
		}

	private:
		void check_piercing_nodes_disjoint() const;
		void check_reachable_sides_are_disjoint() const;
		void check_is_st_cut(nodeid s, nodeid t);
		void check_flow_conservation() const;
		void check_fc_invariants() const;
		void compare_mixed_hyperedges_with_cut() const;
		void check_every_cut_hyperedge_has_source_pin() const;
		void check_free_nodes_are_not_reachable() const;
		void check_border_vertices_are_not_assimilated() const;
	};


	void HLAFlowCutter::check_free_nodes_are_not_reachable() const {
#ifndef NDEBUG
		for (nodeid fn : free_nodes) {
			assert(number_of_mixed_incident_hyperedges[fn] == hg.degree(fn));
			assert(!isSource(fn)); assert(!isTarget(fn));
			assert(!isReachableFromSource(fn)); assert(!isReachableFromTarget(fn));
			assert(isFreeNode(fn));
		}
#endif
	}

	void HLAFlowCutter::check_every_cut_hyperedge_has_source_pin() const {
#ifndef NDEBUG
		for (netid e : source_cut_front) {
			bool has_source_pin = false;
			bool has_nonsource_pin = false;
			for (nodeid pin : hg.pinsOf(e)) {
				has_source_pin |= isSource(pin);
				has_nonsource_pin |= !source_nodes.contains(pin);
			}
			assert(has_source_pin); assert(has_nonsource_pin);
		}
		for (netid e : source_mixed_cut) {
			bool has_source_pin = false;
			bool has_nonsource_pin = false;
			bool has_target_pin = false;
			for (nodeid pin : hg.pinsOf(e)) {
				has_source_pin |= isSource(pin);
				has_nonsource_pin |= !source_nodes.contains(pin);
				has_target_pin |= isTarget(pin);
			}
			assert(has_source_pin); assert(has_nonsource_pin);
		}
#endif
	}

	void HLAFlowCutter::check_piercing_nodes_disjoint() const {
#ifndef NDEBUG
		for (nodeid u : source_piercing_nodes) {
			assert(std::find(target_piercing_nodes.begin(), target_piercing_nodes.end(), u) == target_piercing_nodes.end());
		}
#endif
	}

	void HLAFlowCutter::check_reachable_sides_are_disjoint() const {
#ifndef NDEBUG
		for (nodeid u = 0; u < hg.numNodes(); u++) {
			assert(!(source_nodes.contains(u) && target_nodes.contains(u)));
		}
#endif
	}

	void HLAFlowCutter::compare_mixed_hyperedges_with_cut() const {
#ifndef NDEBUG
		std::vector<netid> cut_from_partition;
		for (netid e = 0; e < hg.numHyperedges(); e++) {
			bool has_source = false;
			bool has_other = false;
			for (nodeid pin : hg.pinsOf(e)) {
				if (source_nodes.contains(pin)) { has_source = true; }
				else { has_other = true; }
			}
			if (has_source && has_other) {
				cut_from_partition.push_back(e);
			}
		}
		std::vector<netid> sort_cut(source_cut_front);
		sort_cut.insert(sort_cut.end(), source_mixed_cut.begin(), source_mixed_cut.end());
		std::sort(sort_cut.begin(), sort_cut.end());
		assert(source_cut_front.size() + source_mixed_cut.size() == cut_from_partition.size());
		for (std::size_t i = 0; i < cut_from_partition.size(); i++) {
			assert(cut_from_partition[i] == sort_cut[i]);
		}
#endif
	}

	void HLAFlowCutter::check_is_st_cut(nodeid s, nodeid t) {
#ifndef NDEBUG
		boost::dynamic_bitset<> is_cut_he(hg.numHyperedges());
		for (netid e : source_cut_front) { is_cut_he.set(e); }
		for (netid e : source_mixed_cut) { is_cut_he.set(e); }
		boost::dynamic_bitset<> node_seen(hg.numNodes()); boost::dynamic_bitset<> he_seen(hg.numHyperedges());
		std::vector<nodeid> queue(hg.numNodes());

		int old_search_direction = source_side_in_current_direction;
		setSearchDirection(0);
		queue[0] = s; node_seen.set(s);
		nodeid qfront = 0; nodeid qend = 1;
		while (qfront < qend) {
			nodeid u = queue[qfront++];
			for (netid e : hg.hyperedgesOf(u)) {
				if (!he_seen[e] && !is_cut_he[e]) {
					he_seen.set(e);
					for (nodeid pin : hg.pinsOf(e)) {
						if (isTarget(pin)) {
							std::cout << pin << std::endl;
							std::cout << isSource(pin) << " " << source_nodes.contains(pin) << " " << hg.degree(pin) << " " << e << std::endl;
						}
						assert(!isTarget(pin));
						assert(!target_nodes.contains(pin));
						if (!node_seen[pin]) {
							queue[qend++] = pin; node_seen.set(pin);
						}
					}
				}
			}
		}
		for (nodeid u = 0; u < hg.numNodes(); u++) {
			if (target_nodes.contains(u)) {
				assert(!node_seen[u]);
			}
		}
		//DO NOT FLIP SEARCH DIRECTION!
		node_seen.reset(); he_seen.reset();
		queue[0] = t; node_seen.set(t);
		qfront = 0; qend = 1;
		while (qfront < qend) {
			nodeid u = queue[qfront++];
			for (netid e : hg.hyperedgesOf(u)) {
				if (!he_seen[e] && !is_cut_he[e]) {
					he_seen.set(e);
					for (nodeid pin : hg.pinsOf(e)) {
						assert(!source_nodes.contains(pin));
						if (!node_seen[pin]) {
							queue[qend++] = pin; node_seen.set(pin);
						}
					}
				}
			}
		}
		for (nodeid u = 0; u < hg.numNodes(); u++) {
			if (source_nodes.contains(u)) {
				assert(!node_seen[u]);
			}
		}
		setSearchDirection(old_search_direction);
#endif
	}

	void HLAFlowCutter::check_flow_conservation() const {
#ifndef NDEBUG
		flow_t src_excess = 0; flow_t tar_excess = 0;
		for (nodeid u = 0; u < hg.numNodes(); u++) {
			flow_t excess = 0;
			for (netid e : hg.hyperedgesOf(u)) {
				if (flow_to[e] == u) excess--;
				if (flow_from[e] == u) excess++;
				assert(!(flow_from[e] == u && flow_to[e] == u));
			}
			if (isSource(u)) { src_excess += excess; }
			if (isTarget(u)) { tar_excess += excess; }
			assert(!(isSource(u) && isTarget(u)));
			if (!isSource(u) && !isTarget(u)) {
				assert(excess == 0 && "flow conservation violated at non-terminal node");
			}
		}
		assert(src_excess == -tar_excess);
#endif
	}

	void HLAFlowCutter::check_fc_invariants() const {
#ifndef NDEBUG
		nodeid src_r=numberOfContractedSourceAssimilatedNodes(), src_a=numberOfContractedSourceAssimilatedNodes();
		nodeid tar_r=numberOfContractedTargetAssimilatedNodes(), tar_a=numberOfContractedTargetAssimilatedNodes();
		for (nodeid u = 0; u < hg.numNodes(); u++) {
			if (isSource(u)) {
				assert(source_nodes.contains(u));
				src_a++;
			}
			if (isTarget(u)) {
				assert(target_nodes.contains(u));
				tar_a++;
			}
			if (isReachableFromSource(u)) {
				src_r++;
			}
			if (isReachableFromTarget(u)) {
				tar_r++;
			}
		}
		assert(src_r >= src_a);
		assert(tar_r >= tar_a);
		assert(numberOfSourceAssimilatedNodes() >= 1);
		assert(numberOfTargetAssimilatedNodes() >= 1);
		assert(src_r == numberOfSourceReachableNodes());
		assert(src_a == numberOfSourceAssimilatedNodes());
		assert(tar_r == numberOfTargetReachableNodes());
		assert(tar_a == numberOfTargetAssimilatedNodes());
		check_flow_conservation();
#endif
	}

	void HLAFlowCutter::check_border_vertices_are_not_assimilated() const {
#ifndef NDEBUG
		for (nodeid p : source_border_vertices) { assert(!isSource(p)); assert(!isTarget(p)); }
#endif
	}

}
