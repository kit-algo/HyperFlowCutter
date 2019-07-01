#pragma once

#include <tlx/unused.hpp>
#include "../extern/patoh_interface.h"
#include "../datastructure/pareto.h"
#include "../util/union_hypergraphs.h"
#include "../util/induced_sub_hypergraph.h"


namespace hyper {

namespace Regrow {
	struct CandidateSet {
		nodeid union_component_to_generate_terminal_for = 0;
		nodeid max_block_size_in_component = 0;
		std::vector<STPair> candidates;
		std::vector<bool> __s_is_random, __t_is_random;
		bool s_is_random(size_t idx) const { return __s_is_random[idx]; }
		bool t_is_random(size_t idx) const { return __t_is_random[idx]; }
	};

	class RegrowBase {
	public:
		RegrowBase(const RegrowBase&) = delete;
		RegrowBase(RegrowBase&&) = delete;
		RegrowBase& operator= (const RegrowBase&) = delete;
		RegrowBase& operator= (RegrowBase&&) = delete;
		virtual ~RegrowBase() = default;

		virtual std::vector<nodeid> regrowBlock(const Hypergraph& union_gap_cc_hgs, std::vector<nodeid>& nodes_in_block, InBlockPredicate& is_in_block, nodeid max_block_size) = 0;

	protected:
		RegrowBase() = default;	//only let subclasses instantiate
	};

	class Bisect : public RegrowBase {
	private:
		//target_weight \in [0,1]
		double max_eps_atboundarybetweenproducing2blocksand1block(const double target_weight) const {
			return std::min(target_weight / (1.0 - target_weight), (1.0 - target_weight)/target_weight);
		}

		std::pair<double, double> min_max_target_weight_range(const double target_weight, const double epsilon) const {
			auto max = target_weight * (1.0+epsilon);
			auto other_target_weight = 1.0 - target_weight;
			auto other_max = other_target_weight * (1.0+epsilon);
			auto min = 1.0 - other_max;
			auto direct_min_formula = target_weight * (1.0+epsilon) - epsilon; tlx::unused(direct_min_formula);
			assert(std::abs( min - direct_min_formula ) < 0.00001);
			return std::make_pair(min, max);
		};

		std::pair<double, nodeid> configure_epsilon_and_max_block_size(const Hypergraph& hg, std::vector<nodeid>& nodes_in_block,
																	   InBlockPredicate& is_in_block, nodeid max_block_size) {
			double min_eps = 0.0;
			nodeid modified_max_block_size = max_block_size;
			nodeid perfect_balance = hg.totalNodeWeight() / 2;

			auto eps_candidates = { 0.03, 0.02, 0.01 };
			for (double eps_candidate : eps_candidates) {
				if (max_block_size > (1.0 + eps_candidate) * perfect_balance) {
					modified_max_block_size = static_cast<nodeid>(std::floor(max_block_size / (1.0 + eps_candidate)));
					return std::make_pair(eps_candidate, modified_max_block_size);
				}
			}

			//even if mbs < n/2
			return std::make_pair(min_eps, modified_max_block_size);
		}
	public:
		std::vector<nodeid> regrowBlock(const Hypergraph& union_gap_cc_hgs, std::vector<nodeid>& nodes_in_block, InBlockPredicate& is_in_block, nodeid max_block_size) override {
			std::vector<nodeid> ejected_nodes;
			auto [refine_again, refine_again_to_hg] = InducedSubHypergraph::extractInducedSubHypergraph(union_gap_cc_hgs, [&](nodeid x) { return is_in_block.inBlock(x); });
			auto [epsilon, modified_max_block_size] = configure_epsilon_and_max_block_size(refine_again, nodes_in_block, is_in_block, max_block_size);
			std::vector<int> one_seed = Random::randomSequence<int>(1);
			Partition refined = PaToHInterface::bipartitionWithMaxSmallerBlocksizeWithPatoh(refine_again, one_seed, modified_max_block_size, epsilon)[0];

			auto [max_refined_id, block_members] = refined.membersOfAllSubsets();
			if (block_members.size() != 2 || max_refined_id != 1) throw std::runtime_error("Only one block from regrow through bisect. This means eps and max_block_size are improperly configured.");
			partitionid larger_block = block_members[0].size() > block_members[1].size() ? 0 : 1;
			nodeid ejected_block = block_members[larger_block].size() <= max_block_size ? 1 - larger_block : larger_block;
			auto rem_size = static_cast<nodeid>(block_members[1 - ejected_block].size());
			std::transform(block_members[ejected_block].begin(),block_members[ejected_block].end(),
						   std::back_inserter(ejected_nodes), [id_map=refine_again_to_hg](nodeid x) { return id_map[x]; });

			if (rem_size > max_block_size) throw std::runtime_error("Remaining block larger than requested. remaining: " + std::to_string(rem_size) + " requested: " + std::to_string(max_block_size));
			return ejected_nodes;
		}
	};

	class EjectRandomFromBorder : public RegrowBase {
	public:
		std::vector<nodeid> regrowBlock(const Hypergraph& union_gap_cc_hgs, std::vector<nodeid>& nodes_in_block, InBlockPredicate& is_in_block, nodeid max_block_size) override {
			if (max_block_size > nodes_in_block.size())	throw std::runtime_error("Regrow through eject. max block size > number of nodes in block. regrow unnecessary.");
			std::vector<nodeid> ejected_nodes;
			if (max_block_size == nodes_in_block.size()) return ejected_nodes;
			auto nEjectedNodes = static_cast<nodeid>(nodes_in_block.size() - max_block_size);
			SparseResetBitvector<nodeid> border_nodes(union_gap_cc_hgs.numNodes());

			{	//find border nodes
				for (netid e = 0; e < union_gap_cc_hgs.numHyperedges(); e++)
					if (std::any_of(union_gap_cc_hgs.pinsOf(e).begin(), union_gap_cc_hgs.pinsOf(e).end(), [&](nodeid x) { return !is_in_block.inBlock(x); }))
						for (nodeid x : union_gap_cc_hgs.pinsOf(e))
							if (is_in_block.inBlock(x) && !border_nodes.contains(x))
								border_nodes.add(x);
			}

			{	//shuffle, in case border nodes suffice. we take full neighborhoods (except for one node) by resizing the queue to appropriate size.
				auto& b = border_nodes.get_contained_ids();
				std::shuffle(b.begin(), b.end(), Random::getRNG());	//in case border nodes already exceed max_block_size, we want a random sample.
			}

			{	//random BFS from border nodes. sample random nodes from front and scan.
				boost::dynamic_bitset<> visited = border_nodes.extractBitset();
				boost::dynamic_bitset<> he_visited(union_gap_cc_hgs.numHyperedges());
				LayeredQueue<nodeid> Q(0);
				Q.inject(border_nodes.extractIds(), union_gap_cc_hgs.numNodes());
				Q.finishNextLayer();
				while (Q.numberOfPushedElements() < nEjectedNodes && !Q.empty()) {
					if (Q.currentLayerEmpty())
						Q.finishNextLayer();
					else {
						//draw random node from current layer and scan
						nodeid pos = Random::randomNumber(Q.qfront, Q.layerend - 1);
						nodeid u = Q.swapFrontToPositionAndPop(pos);
						for (netid e : union_gap_cc_hgs.hyperedgesOf(u))
							if (!he_visited[e]) {
								he_visited.set(e);
								for (nodeid p : union_gap_cc_hgs.pinsOf(e))
									if (!visited[p] && is_in_block.inBlock(p)) {
										Q.push(p);
										visited.set(p);
									}
							}
					}
				}
				auto nRequiredNodesFromNonBorderCCs = nEjectedNodes - std::min(nEjectedNodes, Q.numberOfPushedElements());
				bool need_to_sample_from_ccs_without_border = Q.empty() && nRequiredNodesFromNonBorderCCs > 0;
				ejected_nodes = Q.extract();
				ejected_nodes.erase(ejected_nodes.begin() + nEjectedNodes, ejected_nodes.end());
				if (need_to_sample_from_ccs_without_border) {
					//fill up with random nodes. this means they are taken from a cc without a border. so we might as well go ham with randomness.
					//LOG << "border did not suffice. sample additionally";
					std::vector<nodeid> candidates;
					for (nodeid x : nodes_in_block)
						if (!visited[x])
							candidates.push_back(x);
					std::sample(candidates.begin(), candidates.end(), std::back_inserter(ejected_nodes), nRequiredNodesFromNonBorderCCs, Random::getRNG());
				}
			}

			assert(ejected_nodes.size() == nEjectedNodes);
			return ejected_nodes;
		}
	};

class FBT {
public:
	struct MultiCutterOptions {
		bool use_external_partitions_in_result_front = false;
		RegrowBlockStrategy regrow_block_strategy = RegrowBlockStrategy::BisectStrat;
//		RegrowBlockStrategy regrow_block_strategy = RegrowBlockStrategy::RegrowFromRandomNode;
		FBTPartitionStrategy partition_strategy = FBTPartitionStrategy::TwoWay;
//		PartitionStrategy partition_strategy = PartitionStrategy::ThreeWay;
	};

	struct DisconnectedMultiCutterOptions {
		nodeid first_cc_with_gap = 0;
		nodeid cc_sizes_at_first_gap = 0;
	};

	DisconnectedMultiCutterOptions dmco;
	MultiCutterOptions mco;
	std::unique_ptr<Regrow::RegrowBase> regrow_strategy;

private:
	static constexpr bool debug = true;
	ConnectedComponents& cc;
	CandidateSet cs;
	nodeid regrow_threshold = 0;
	nodeid globalNumNodes;
	bool initialized = false;
	bool use_finish_balance_partitioning = true;
	bool s_is_random(nodeid union_component) const { return cs.s_is_random(union_component); }
	bool t_is_random(nodeid union_component) const { return cs.t_is_random(union_component); }

	nodeid union_component_to_generate_terminal_for = 0;
	std::vector<nodeid> prefix_sum_for_local2union;

	flow_t perfectly_balanced_cut_of_external_partitioner = MAX_FLOW;

public:
	FBT(nodeid numNodes, ConnectedComponents& cc) : cc(cc), globalNumNodes(numNodes) {

	};

	void selectComponentForTerminalGeneration(nodeid c) {
		if (c >= cc.numComponents())
			throw std::runtime_error("Unknown component c=" + std::to_string(c) + ".");
		union_component_to_generate_terminal_for = c - dmco.first_cc_with_gap;
	}

	void generateFBTWithPatoh(const Hypergraph& global_hg, const Hypergraph& union_gap_cc_hgs, ParetoFront& result_front, STOptions& stOptions, State& state) {
		if (!stOptions.useFBT()) { return; }
		init(stOptions);
		partitionid unfixable_block = invalid_block;
		auto[intersection_partition, partitions, unfixable_nodes] = generatePartitions(global_hg, union_gap_cc_hgs, state, result_front, stOptions, unfixable_block);
		LOG << "partitioning done";
		LOG << "partition sizes";
		for (auto& x : intersection_partition.getSetSizes())
			LOG << "set size" << x;
		if (intersection_partition.maxIdIfNotModified() < 1) {
			throw std::runtime_error("Less than two blocks in FBP.");
		}
		regrowTooLargeBocks(union_gap_cc_hgs, intersection_partition, partitions, unfixable_nodes, unfixable_block);
		LOG << "regrow too large blocks done";
		setTwoLargestBlocksAsCandidates(intersection_partition, unfixable_block);
		LOG << "candidate selection done";
		for (nodeid union_c = 0; union_c < cc.numComponents() - dmco.first_cc_with_gap; union_c++)
			LOG << cs.candidates[union_c].s.size() << cs.candidates[union_c].t.size() << " size of s and t in union component" << union_c;
		stOptions.perfectly_balanced_cut_with_external_partitioner = perfectly_balanced_cut_of_external_partitioner;
	}

	void enforceSolvability(STOptions& stOptions, const Hypergraph& global_hg, const Hypergraph& union_gap_cc_hgs) {
		if (!stOptions.useFBT()) return;
		Partition dummy_partition(union_gap_cc_hgs.numNodes());
		InBlockPredicate dummy_in_block_predicate(dummy_partition, 1);
		for (nodeid c = dmco.first_cc_with_gap, union_c = 0; c < cc.numComponents(); c++, union_c++) {
			auto& st = cs.candidates[union_c];
			nodeid max_block_size = cc.componentSize(c) / 2;	//maybe be more aggressive here?
			if ( (st.s.size() == cc.componentSize(c) && !s_is_random(union_c)) || (st.t.size() == cc.componentSize(c) && !t_is_random(union_c)) ) {
				cs.__s_is_random[union_c] = true;
				cs.__t_is_random[union_c] = true;
				if (st.s.size() != cc.componentSize(c)) { st.s.resize(cc.componentSize(c)); std::iota(st.s.begin(), st.s.end(), 0); }
				if (st.t.size() != cc.componentSize(c)) { st.t.resize(cc.componentSize(c)); std::iota(st.t.begin(), st.t.end(), 0); }
			}
			else {
				if (!cs.s_is_random(union_c) && st.s.size() > max_block_size) {
					regrowBlockForSolvability(union_gap_cc_hgs, union_c, st.s, dummy_partition, dummy_in_block_predicate, max_block_size);
				}
				if (!cs.t_is_random(union_c) && st.t.size() > max_block_size) {
					regrowBlockForSolvability(union_gap_cc_hgs, union_c, st.t, dummy_partition, dummy_in_block_predicate, max_block_size);
				}
			}
		}
	}

	//this creates copies!
	std::vector<STPair> generateTerminalPairs(const STOptions& stOptions, std::size_t generatorID) {
		if (!stOptions.useFBT()) { return std::vector<STPair>(0); }
		LOG << "generate terminals for " << V(union_component_to_generate_terminal_for) << V(dmco.first_cc_with_gap) << V(unioncc2globalcc(union_component_to_generate_terminal_for));
		sanity_check();
		STPair& st_ref = cs.candidates[union_component_to_generate_terminal_for];
		nodeid union_c = union_component_to_generate_terminal_for;
		nodeid c = unioncc2globalcc(union_c);
		if ( (st_ref.s.size() == cc.componentSize(c) && !s_is_random(union_c)) || (st_ref.t.size() == cc.componentSize(c) && !t_is_random(union_c)) ) {
			LOG << "skip component";
			return std::vector<STPair>(0);	//no cuts necessary
		}
		std::vector<STPair> result(stOptions.numRepetitionsPerFinishBalanceTerminal);
		bool t_random = t_is_random(union_c), s_random = s_is_random(union_c);
		if (t_random && s_random) {		//make the assertion much more elegant please. std::equal changes signatures between C++17 and C++20 :(
#ifndef NDEBUG
			//signature for std::equal changes between current C++ standards --> not pretty
			std::sort(st_ref.s.begin(), st_ref.s.end()); std::sort(st_ref.t.begin(), st_ref.t.end());	//should be sorted already, I think
			assert(st_ref.s.size() == st_ref.t.size());
			for (std::size_t i = 0; i < st_ref.s.size(); i++) assert(st_ref.s[i] == st_ref.t[i]);
#endif
			LOG << "s and t random";
			for (auto& x : result) {
				auto sample = Random::sampleTwoDisjoint(st_ref.s);
				x.s = { sample.first }; x.t = { sample.second };
			}
		}
		else {
			for (auto& x : result) {
				if (s_random) x.s = { Random::sampleOne(st_ref.s) };
				else x.s = st_ref.s;	//make a copy. if we constructed result vector with copies of st_ref, they would get copied anyways.
				if (t_random) x.t = { Random::sampleOne(st_ref.t) };
				else x.t = st_ref.t;	//make a copy.

				if (s_random) LOG << "s random";
				else LOG << "s not random";
				if (t_random) LOG << "t random";
				else LOG << "t not random";
				LOG << V(x.s.size()) << V(x.t.size());
			}
		}
		for (auto& x : result) {
			x.origin = STMethod::FBT;
			x.generatorID = generatorID;
		}
		return result;
	}

private:
	bool connected_hypergraph() const { return cc.numComponents() <= 1; }

	inline nodeid local2union(nodeid u, nodeid union_c) { return prefix_sum_for_local2union[union_c] + u; }

	inline nodeid local_nodeid_2union_cc_id(nodeid u) {
		assert (u < prefix_sum_for_local2union.back());
		auto iter = std::lower_bound(prefix_sum_for_local2union.begin(), prefix_sum_for_local2union.end(), u+1) - 1;	//psfl2u[0] = 0 --> u+1 ensures iter points at least to psfl2u[1]
		return static_cast<nodeid>(std::distance(prefix_sum_for_local2union.begin(), iter));
	}

	inline nodeid union2local(nodeid u, nodeid union_c) {
		return u - prefix_sum_for_local2union[union_c];
	}

	inline nodeid union2local(nodeid u) {
		nodeid union_c = local_nodeid_2union_cc_id(u);
		return union2local(u, union_c);
	}

	inline nodeid unioncc2globalcc(nodeid union_c) { return union_c + dmco.first_cc_with_gap; }
	inline nodeid glocalcc2unioncc(nodeid c) { return c - dmco.first_cc_with_gap; }

	void init(STOptions& stOptions) {
		initialized = true;
		mco.regrow_block_strategy = stOptions.fbt_regrow_strat;
		mco.partition_strategy = stOptions.fbt_partition_strat;

		regrow_threshold = static_cast<nodeid>(std::floor(stOptions.fbt_fractional_blocksize * globalNumNodes));
		if (stOptions.fbt_fractional_blocksize == 0.5)
			regrow_threshold = static_cast<nodeid>(globalNumNodes/2);	//do it exact for the floating point imprecision haters


		auto nComp = static_cast<nodeid>(cc.numComponents() - dmco.first_cc_with_gap);
		cs.candidates = std::vector<STPair>(nComp);
		cs.__s_is_random = std::vector<bool>(nComp, false);
		cs.__t_is_random = std::vector<bool>(nComp, false);
		prefix_sum_for_local2union = std::vector<nodeid>(nComp + 1, 0);
		std::partial_sum(cc.ccsizes.begin() + dmco.first_cc_with_gap, cc.ccsizes.end(), prefix_sum_for_local2union.begin() + 1);

		switch (mco.regrow_block_strategy) {
			case RegrowBlockStrategy::BisectStrat :
				regrow_strategy = std::make_unique<Regrow::Bisect>();
				break;
			case RegrowBlockStrategy ::EjectRandomFromBorderStrat:
				regrow_strategy = std::make_unique<Regrow::EjectRandomFromBorder>();
				break;
			default:
				std::cout << "Unknown regrow block partition_strategy" << std::endl;
		}
	}

	void sanity_check() {
		if (!initialized) throw std::runtime_error("FinishBalancePartitioning not initialized");
		if (!use_finish_balance_partitioning) throw std::runtime_error("FinishBalancePartitioning not intended for use.");
	}

	void if_activated_then_write_external_partitions_to_pareto_front(const Hypergraph &global_hg,
																	 const Hypergraph &cc_hg,
																	 ParetoFront &result_front,
																	 std::vector<Partition> &partitions,
																	 State& state)
	{
		nodeid global_perf_balance = global_hg.totalNodeWeight() / 2;
		for (auto& p : partitions) {
			nodeid a0 = p.getSetSize(0);
			nodeid b0 = p.getSetSize(1);
			auto [mi,ma] = std::minmax(a0,b0);
			flow_t cut = p.computeCut(cc_hg);
			nodeid reach_with_smaller_block = std::min(global_perf_balance, mi + dmco.cc_sizes_at_first_gap);
			nodeid best_balanced_smaller_block = reach_with_smaller_block;

			std::cout << "PaToH cut: " << cut << " mi " << mi << " ma " << ma << std::endl;

			if (mco.use_external_partitions_in_result_front && cut < result_front.cuts[reach_with_smaller_block].cut) {
				result_front.cuts[reach_with_smaller_block] =
						CutInfo::getCutFromExternalPartitioner(cut, reach_with_smaller_block,
															   global_hg.totalNodeWeight() - reach_with_smaller_block, duration(0.0));
			}
			if (reach_with_smaller_block == global_perf_balance)
				perfectly_balanced_cut_of_external_partitioner = std::min(perfectly_balanced_cut_of_external_partitioner, cut);
			if (ma <= global_perf_balance) {
				nodeid reach_with_larger_block = std::min(global_perf_balance, ma + dmco.cc_sizes_at_first_gap);
				best_balanced_smaller_block = std::max(best_balanced_smaller_block, reach_with_larger_block);
				if (mco.use_external_partitions_in_result_front && cut < result_front.cuts[reach_with_larger_block].cut) {
					result_front.cuts[reach_with_smaller_block] =
							CutInfo::getCutFromExternalPartitioner(cut, reach_with_larger_block,
																   global_hg.totalNodeWeight() - reach_with_larger_block, duration(0.0));

				}
				if (reach_with_larger_block)
					perfectly_balanced_cut_of_external_partitioner = std::min(perfectly_balanced_cut_of_external_partitioner, cut);
			}

			state.external_partitioner_results.push_back( {cut, best_balanced_smaller_block} );
		}
	}

	std::tuple<Partition, std::vector<Partition>, std::vector<bool>> generatePartitions(const Hypergraph& global_hg, const Hypergraph& union_gap_cc_hgs, State& state,
																						ParetoFront& result_front, STOptions& stOptions, partitionid& unfixable_block) {
		// not such a great idea, maybe. maybe do half or 3/4 of the non-gap nodes.
		nodeid max_smaller_blocksize = global_hg.totalNodeWeight()/2 - dmco.cc_sizes_at_first_gap; 	//aim for shuffling all non-gap components to the same side to get the smallest possible cut with PaToH. Hope for epsilon to sort of mitigate the rest.
		//if cc_sizes_at_first_gap = 0. does that mean we should set epsilon = 0? rather not. give epsilon a base value, which we then increase. the base value is a user-parameter.
		std::vector<int> seeds = Random::randomSequence<int>(stOptions.numExternalPartitionerCallsPerFinishBalanceTerminal);
		std::vector<Partition> partitions;
		std::vector<bool> unfixable_nodes(union_gap_cc_hgs.numNodes(), false);
		if (mco.partition_strategy == FBTPartitionStrategy::TwoWay) {
			LOG << V(stOptions.fbt_epsilon);
			partitions = PaToHInterface::bipartitionWithMaxSmallerBlocksizeWithPatoh(union_gap_cc_hgs, seeds, max_smaller_blocksize, stOptions.fbt_epsilon, stOptions.patoh_preset);
			if_activated_then_write_external_partitions_to_pareto_front(global_hg, union_gap_cc_hgs, result_front, partitions, state);
		}
		else if (mco.partition_strategy == FBTPartitionStrategy::ThreeWay) {	//This setting makes little sense when using multiple external partitioner calls for confidence boosting. still keep it :)
			partitions = PaToHInterface::threeWayPartitionWithTwoEqualSizedBlocksWithPatoh(global_hg, seeds);
			//make sure the third (smallest block) does not get used for preassignment?
			for (auto& p : partitions) {
				auto m = p.membersAndIdsOfAllSubsets();
				auto third_block_it = std::min_element(m.begin(), m.end(), [](const auto& lhs, const auto& rhs) { return lhs.second.size() < rhs.second.size(); });
				for (nodeid u : third_block_it->second) { unfixable_nodes[u] = true; }
			}
		}
		else {
			throw std::runtime_error("Unknown FBT partition strategy.");	//dead code
		}

		Partition intersection_partition = PartitionIntersection::bulkIntersect(partitions);
		unfixable_block = intersection_partition.maxIdIfNotModified() + 1;
		LOG << V(unfixable_block) << V(intersection_partition.maxId());
		if (unfixable_block < 2)
			throw std::runtime_error("Unfixable block id too small because too few blocks.");
		for (nodeid u = 0; u < union_gap_cc_hgs.numNodes(); u++)
			if (unfixable_nodes[u])
				intersection_partition[u] = unfixable_block;
		return std::make_tuple(intersection_partition, partitions, unfixable_nodes);
	};

	void regrowTooLargeBlock(const Hypergraph& union_gap_cc_hgs, Partition& intersection_partition, std::vector<Partition>& partitions,
							 std::vector<bool>& unfixable_nodes, partitionid& unfixable_block, std::vector<nodeid>& nodes_in_block, partitionid p, nodeid max_block_size) {
		InBlockPredicate ibp(intersection_partition, p);
		auto ejected_candidates = regrow_strategy->regrowBlock(union_gap_cc_hgs, nodes_in_block, ibp, max_block_size);
		LOG << "regrow block p=" << p <<  V(max_block_size) << V(nodes_in_block.size()) << V(ejected_candidates.size());
		for (nodeid x : ejected_candidates) {
			intersection_partition[x] = unfixable_block;
			unfixable_nodes[x] = true;
		}
	}

	void regrowTooLargeBocks(const Hypergraph& union_gap_cc_hgs, Partition& intersection_partition, std::vector<Partition>& partitions, std::vector<bool>& unfixable_nodes, partitionid& unfixable_block) {
		for (auto&[p, m] : intersection_partition.membersAndIdsOfAllSubsets())
			if (m.size() > regrow_threshold) regrowTooLargeBlock(union_gap_cc_hgs, intersection_partition, partitions, unfixable_nodes, unfixable_block, m, p, regrow_threshold);
	}

	void setTwoLargestBlocksAsCandidates(Partition& intersection_partition, partitionid unfixable_block) {
		partitionid max_id = intersection_partition.maxId();
		for (nodeid c = dmco.first_cc_with_gap, union_c = 0; c < cc.numComponents(); c++, union_c++) {
			LOG << "select for " << V(c) << V(union_c);
			auto& st = cs.candidates[union_c];
			std::vector< std::pair<partitionid, std::vector<nodeid> >> overlap(max_id + 1);
			partitionid i = 0; for (auto& x : overlap) { x.first = i++; }
			for (nodeid u_local = 0; u_local < cc.componentSize(c); u_local++)
				overlap[ intersection_partition[ local2union(u_local, union_c) ] ].second.push_back(u_local);
			std::sort(overlap.begin(), overlap.end(), [&](const auto& lhs, const auto& rhs) { return lhs.second.size() < rhs.second.size(); });
			assert(overlap.size() >= 2);

			auto select_block = [&](std::vector<nodeid>& cand, std::vector<bool>& is_random) {
				if (overlap.back().first == unfixable_block)
					overlap.pop_back();
				if (overlap.back().second.size() <= 1) {
					is_random[union_c] = true;
					cand.resize(cc.componentSize(c)); std::iota(cand.begin(), cand.end(), 0);
				}
				else
					cand = std::move(overlap.back().second);
				overlap.pop_back();
			};
			select_block(st.s, cs.__s_is_random);
			select_block(st.t, cs.__t_is_random);

			if (t_is_random(union_c) && !s_is_random(union_c)) {
				//this is unnecessary if PaToH always produced at least two blocks. However PaToH may put everything in one block, when we have many vertices in non_gap ccs
				//eject s from t's candidates
				std::vector<nodeid> t_cands;
				std::set_difference(st.t.begin(), st.t.end(), st.s.begin(), st.s.end(), std::back_inserter(t_cands));
				std::swap(t_cands, st.t);

				LOG << "ejected candidates from t." << st.t.size() << " is the new number of candidates";
			}

#ifndef NDEBUG
			std::sort(st.s.begin(), st.s.end()); std::sort(st.t.begin(), st.t.end());
			std::vector<nodeid> s_t_intersection;
			std::set_intersection(st.s.begin(), st.s.end(), st.t.begin(), st.t.end(), std::back_inserter(s_t_intersection));
			assert(cs.__s_is_random[union_c] || cs.__t_is_random[union_c] || s_t_intersection.empty());
#endif

		}
	}

	void regrowBlockForSolvability(const Hypergraph& union_gap_cc_hgs, nodeid union_c, std::vector<nodeid>& nodes_in_block, Partition& dummy_partition, InBlockPredicate& dummy_in_block_predicate, nodeid max_block_size) {
		dummy_in_block_predicate.blockid++;
		nodeid c = unioncc2globalcc(union_c);
		std::vector<nodeid> nodes_in_block_union_ids;
		for (nodeid x : nodes_in_block) {
			nodeid union_x = local2union(x, union_c);
			dummy_partition[union_x] = dummy_in_block_predicate.blockid;	//enable in predicate
			nodes_in_block_union_ids.push_back(union_x);
		}
		auto ejected_nodes = regrow_strategy->regrowBlock(union_gap_cc_hgs, nodes_in_block_union_ids, dummy_in_block_predicate, max_block_size);
		std::vector<bool> is_ejected_local(cc.componentSize(c), false);
		for (nodeid union_x : ejected_nodes) is_ejected_local[ union2local(union_x, union_c) ] = true;
		auto remove_it = std::remove_if(nodes_in_block.begin(), nodes_in_block.end(), [&](nodeid local_x) { return is_ejected_local[local_x]; });
		nodes_in_block.erase(remove_it, nodes_in_block.end());
	}

};

}//namespace Regrow
}//namespace hyper
