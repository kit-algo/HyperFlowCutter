#pragma once

#include "../datastructure/id_queue.hpp"
#include "../datastructure/hypergraph.h"
#include "hla_flowcutter.h"
#include "hypergraph_search_algorithms.h"
#include "../datastructure/queue.h"
#include "connectivity_check.h"
#include "hla_piercing_heuristics/patoh_ensemble_partitioning.h"
#include "terminals.h"
#include "../datastructure/pareto.h"
#include "hla_vertex_disjoint_edmonds_karp.h"
#include "hla_dinic.h"


namespace hyper {
	class MultiCutter {
	public:
		template<typename PierceScorer, typename FTreatCut>
		static void runInterleaved(const Hypergraph &hg, const std::vector<STPair> &st_pairs, SearchAlgorithm& search_algorithm,
								   PierceScorer scorer, OverlapWithPartitionBlocks& opb, FTreatCut treatCut, nodeid max_block_size, flow_t external_cut_bound) {


			bool copy_initial_cut = std::all_of(st_pairs.begin(), st_pairs.end(), [](const STPair& st) { return st.origin == STMethod::FBT && st.s_is_coarse && st.t_is_coarse; });
			//in this case we can guarantee the initial cut is equal among all cutters. so a copy should be faster than computing from scratch

			//std::vector<HLAFlowCutter> cutters(st_pairs.size(), HLAFlowCutter(hg, search_algorithm, opb));
			std::cout << "external cut bound= " << external_cut_bound << std::endl;
			std::vector<HLAFlowCutter> cutters;
			cutters.reserve(st_pairs.size());

			flow_t best_desiredly_balanced_cutsize = external_cut_bound;
			MinIDQueue active_cutters(static_cast<nodeid>(st_pairs.size()));
			uint16_t i = 0;

			auto create_cutter = [&](const STPair& st) {
				cutters.emplace_back( hg, search_algorithm, opb );
				cutters.back().desired_smaller_blocksize = max_block_size;
				cutters.back().initialCut(st);
				if (cutters.back().isDesiredlyBalanced()) {
					best_desiredly_balanced_cutsize = std::min(best_desiredly_balanced_cutsize, cutters.back().cutsize);
				}
				treatCut(cutters.back(), i);
				if (cutters.back().cutsize < best_desiredly_balanced_cutsize) { active_cutters.push({i, cutters.back().cutsize}); }
				assert(i < st_pairs.size());
				i++;
			};

			if (copy_initial_cut) {
				const STPair& st = st_pairs.front();
				create_cutter(st);
				const auto& first_cutter = cutters.front();
				for (std::size_t j = 1; j < st_pairs.size(); j++) {
					cutters.push_back( first_cutter );
					if (cutters.back().cutsize < external_cut_bound) { active_cutters.push({i, cutters.back().cutsize}); }
					assert(i < st_pairs.size());
					i++;
				}
			}
			else {
				for (auto& st : st_pairs) {
					create_cutter(st);
				}
			}



			while (!active_cutters.empty() && active_cutters.peek().key <= best_desiredly_balanced_cutsize) {
				nodeid mincutter_id = active_cutters.peek().id;
				auto &c = cutters[mincutter_id];
				if (!c.advance(scorer) || c.cutsize >= best_desiredly_balanced_cutsize) {
					active_cutters.pop();
					if (c.cutsize == best_desiredly_balanced_cutsize)
						treatCut(c, mincutter_id);
				}
				else {
					treatCut(c, mincutter_id);
					if (c.isDesiredlyBalanced()) {
						best_desiredly_balanced_cutsize = std::min(best_desiredly_balanced_cutsize, c.cutsize);
					}
					active_cutters.increase_key({mincutter_id, c.cutsize});
				}
			}
		}

		template<typename PierceScorer, typename FTreatCut>
		static void runConsecutive(const Hypergraph& hg, const std::vector<STPair>& st_pairs, SearchAlgorithm& search_algorithm,
								   PierceScorer& scorer, OverlapWithPartitionBlocks& opb, FTreatCut treatCut, nodeid max_block_size, flow_t external_cut_bound) {
			uint16_t i = 0;
			for (auto& st : st_pairs) {
				HLAFlowCutter cutter(hg, search_algorithm, opb);
				cutter.desired_smaller_blocksize = max_block_size;
				cutter.initialCut(st);
				treatCut(cutter, i);

				while (!cutter.isDesiredlyBalanced() && cutter.advance(scorer)) {
					treatCut(cutter, i);
				}
				i++;
			}
		}

		template<typename PierceScorer, typename FTreatCut>
		static void runHybrid(const Hypergraph& hg, const std::vector<STPair>& st_pairs, SearchAlgorithm& search_algorithm,
							  PierceScorer scorer, OverlapWithPartitionBlocks& opb, FTreatCut treatCut, nodeid max_block_size, flow_t external_cut_bound) {
			uint32_t nFlowIncreasingIterations = 10;
			flow_t flowGap = 150;
			//std::vector<HLAFlowCutter> cutters(st_pairs.size(), HLAFlowCutter(hg, search_algorithm, opb));
			std::vector<HLAFlowCutter> cutters;
			cutters.reserve(st_pairs.size());
			for (std::size_t i = 0; i < st_pairs.size(); i++) { cutters.emplace_back( hg, search_algorithm, opb ); }	//lack of copy constructor

			flow_t best_desiredly_balanced_cutsize = std::numeric_limits<flow_t>::max();
			MinIDQueue active_cutters(static_cast<nodeid>(st_pairs.size()));

			uint16_t i = 0;
			for (auto& st : st_pairs) {
				cutters[i].desired_smaller_blocksize = max_block_size;
				cutters[i].initialCut(st);
				if (cutters[i].isDesiredlyBalanced()) {
					best_desiredly_balanced_cutsize = std::min(best_desiredly_balanced_cutsize, cutters[i].cutsize);
				}
				treatCut(cutters[i], i);
				if (cutters[i].cutsize < external_cut_bound) { active_cutters.push({i, cutters[i].cutsize}); }
				i++;
			}

			while (!active_cutters.empty() && active_cutters.peek().key <= best_desiredly_balanced_cutsize) {
				nodeid mincutter_id = active_cutters.peek().id;
				auto &c = cutters[mincutter_id];
				flow_t start_flow = c.cutsize;
				flow_t last_flow = start_flow;
				uint32_t iterationsInWhichFlowIncreased = 0;
				bool cutter_cannot_advance = false;
				//take whichever criterion lasts longer
				while (c.cutsize <= start_flow + flowGap || iterationsInWhichFlowIncreased < nFlowIncreasingIterations) {
					if (!c.advance(scorer)) {
						cutter_cannot_advance = true; break;
					}
					else {
						treatCut(cutters[mincutter_id], mincutter_id);
						if (c.isDesiredlyBalanced()) {
							best_desiredly_balanced_cutsize = std::min(best_desiredly_balanced_cutsize, c.cutsize);
							cutter_cannot_advance = true; break;
						}
						if (c.cutsize > last_flow) {
							last_flow = c.cutsize;
							iterationsInWhichFlowIncreased++;
						}
					}
				}
				if (cutter_cannot_advance || c.cutsize >= external_cut_bound) {
					active_cutters.pop();
				}
				else {
					active_cutters.increase_key({mincutter_id, c.cutsize});
				}
			}
		}

		template<class Piercer, typename ExecuteAfterSTOption>
		static void run(const Hypergraph& hg, ParetoFront& front, State& state, EnsemblePartitioning& ep, std::vector<STOptions>& stoptions, nodeid mbs, ExecuteAfterSTOption easto) {
			ConnectedComponents cc = ConnectedComponents::createCCOfConnectedHG(hg.totalNodeWeight());
			flow_t external_cut_bound = MAX_FLOW;
			for (auto& st : stoptions) {
				external_cut_bound = std::min(external_cut_bound, front.getCut(mbs).second.cut);

				auto t = time_now();
				Regrow::FBT fbp(hg.totalNodeWeight(), cc);
				fbp.generateFBTWithPatoh(hg, hg, front, st, state);
				for (auto& x : state.external_partitioner_results) {
					if (x.size_of_smaller_block >= mbs) {
						external_cut_bound = std::min(external_cut_bound, x.cut);
					}
				}
				auto t_fbt = time_now();
				state.running_time_gen_fbt += duration(t_fbt - t);

				t = t_fbt;
				run<Piercer>(hg, front, ep, fbp, state, st, mbs, external_cut_bound);
				state.running_time_multi_cutter += duration(time_now() - t);
				easto(st, hg, front, state);
			}
		}

		//default way to call HyperFlowCutter. if you want extra output on your ParetoFront...call the other function with your custom treatCut lambda
		template<class Piercer>
		static void run(const Hypergraph& hg, ParetoFront& front, EnsemblePartitioning& ep, Regrow::FBT& fbp, State& state, STOptions& stoptions, nodeid mbs, flow_t external_cut_bound) {
			external_cut_bound = std::min(front.getCut(mbs).second.cut, external_cut_bound);	//this also compares the global cut bound with the cc-local one. fortunately the cc-local one is always a lower bound on the global one, so this is fine!

			STChooser stpairfinder(hg);
			auto [st_pairs, furtherAway, farAway, random] = stpairfinder.getSTPairs(hg, state, stoptions, ep, fbp);
			tlx::unused(furtherAway, farAway, random);
			if (st_pairs.empty())
				return;
			//set up pareto front datastructures
			std::vector<ParetoWrites> pw(st_pairs.size(), ParetoWrites(front));
			auto t = time_now();
			auto t_begin_step = t;
			auto treatCut = [&](HLAFlowCutter& flc, int cutter_id) {
				pw[cutter_id].writeAllAvailableCuts(flc, cutter_id, time_now() - t);
				if (state.output_detail >= 200) {
					auto t_end_step = time_now();
					std::cout << cutter_id << ": " << flc.state() << " time: " << inMilliseconds(duration(t_end_step - t_begin_step)).count() << "ms" << "\n";
					t_begin_step = t_end_step;
				}
			};

			std::size_t genID = st_pairs.front().generatorID;
			bool use_coarse_hypergraph = std::all_of(st_pairs.begin(), st_pairs.end(), [genID](const auto& st) { return st.origin == STMethod::FBT && st.generatorID == genID; } );
			use_coarse_hypergraph |= (st_pairs.size() == 1 && st_pairs[0].origin == STMethod::Ensemble);
			if (use_coarse_hypergraph) {
				Hypergraph coarse_hg = contractTerminalPairs(hg, st_pairs);
				run<Piercer>(coarse_hg, st_pairs, state, ep, treatCut, external_cut_bound, mbs);
			}
			else {
				run<Piercer>(hg, st_pairs, state, ep, treatCut, external_cut_bound, mbs);
			}

			std::cout << std::flush;
		}

		template<class Piercer, typename FTreatCut>
		static void run(const Hypergraph& hg, std::vector<STPair>& st_pairs, State& state,
						EnsemblePartitioning& ep, FTreatCut treatCut, flow_t external_cut_bound, nodeid mbs) {

			std::unique_ptr<SearchAlgorithm> search_algorithm = instantiateSearchAlgorithm(hg, state);

			Piercer pp(hg.numHyperedges(), state, ep);
			OverlapWithPartitionBlocks opb(ep.intersectionPartition, state.piercerOptions.use_ensemble);

			if (state.multiCutterExecution == MultiCutterExecution::Interleaved) {
				runInterleaved(hg, st_pairs, *search_algorithm, pp, opb, treatCut, mbs, external_cut_bound);
			}
			else if (state.multiCutterExecution == MultiCutterExecution::Consecutive) {
				runConsecutive(hg, st_pairs, *search_algorithm, pp, opb, treatCut, mbs, external_cut_bound);
			}
			else if (state.multiCutterExecution == MultiCutterExecution::Hybrid) {
				runHybrid(hg, st_pairs, *search_algorithm, pp, opb, treatCut, mbs, external_cut_bound);
			}
		}


	private:
		static std::unique_ptr<SearchAlgorithm> instantiateSearchAlgorithm(const Hypergraph &hg, State &state) {
			if (state.type_of_flow_algorithm == TypeOfFlowAlgorithm::Dinic) return std::make_unique<HLADinic>(hg, state.build_datastructures_during_grow_reachable);
			else if (state.type_of_flow_algorithm == TypeOfFlowAlgorithm::VertexDisjointEdmondsKarp) return std::make_unique<VertexDisjointHLABFS>(hg, state.build_datastructures_during_grow_reachable);
			else throw std::runtime_error("Unkown flow algorithm");
		}


		static Hypergraph contractTerminalPairs(const Hypergraph& fine_hg, std::vector<STPair>& terminals) {
			if (terminals.empty()) throw std::runtime_error("Trying to contract empty set of terminals");
			if (std::any_of(terminals.begin(), terminals.end(), [](const STPair& st) { return st.origin != STMethod::FBT; })
				&& !(terminals.size() == 1 && terminals[0].origin == STMethod::Ensemble)) {
				throw std::runtime_error("Contraction is only supported for FBT terminals or single ensemble terminals.");
			}
			std::vector<std::vector<nodeid>> coarse_nodes;
			bool s_is_coarse = terminals.front().s.size() > 1;
			if (s_is_coarse)
				coarse_nodes.emplace_back(std::move(terminals.front().s));
			bool t_is_coarse = terminals.front().t.size() > 1;
			if (t_is_coarse)
				coarse_nodes.emplace_back(terminals.front().t);
			for (auto& st : terminals) {
				st.s = {0}; st.t = {1};
				if (s_is_coarse)
					st.s = {0};
				if (t_is_coarse)
					st.t = { s_is_coarse ? nodeid(1) : nodeid(0) };
				st.s_is_coarse = s_is_coarse;
				st.t_is_coarse = t_is_coarse;
			}
			return Hypergraph(fine_hg, coarse_nodes);
		}

	};

}
