#pragma once

#include "multicutter.h"
#include "connectivity_check.h"

namespace hyper {



class ExtractConnectedComponents {
public :
	//subhypergraphs and global2local nodeid mapping
	static std::pair<std::vector<Hypergraph>, std::vector<nodeid>> extractConnectedComponentHypergraphs(const Hypergraph& hg, ConnectedComponents& cc) {
		std::vector<nodeid> global2local(hg.numNodes());
		std::vector<nodeid> cc_current_nodeid(cc.numComponents(), 0);
		for (nodeid u = 0; u < hg.numNodes(); u++) { global2local[u] = cc_current_nodeid[cc.node_component[u]]++; }

		std::vector<std::vector<nodeid>> hyperedge_sizes(cc.numComponents());
		std::vector<std::vector<nodeid>> pins(cc.numComponents());
		for (netid e = 0; e < hg.numHyperedges(); e++) {
			nodeid c = cc.hyperedge_component[e];
			hyperedge_sizes[c].push_back(hg.pinCount(e));
			for (nodeid pin : hg.pinsOf(e)) {
				pins[c].push_back(global2local[pin]);
			}
		}

		std::vector<Hypergraph> hgs;
		for (nodeid c = 0; c < cc.numComponents(); c++) {	//never delete large hyperedges here. if the cc gets disconnected again, we have to repeat recursively and distinguish all those in the DP. too annoying.
			hgs.emplace_back(cc.ccsizes[c], hyperedge_sizes[c], std::move(pins[c]));
		}
		return std::make_pair(hgs, global2local);
	}
};


class DisconnectedMultiCutter {
public:
	static constexpr bool debug = false;

	//sort ccs ascending by size
	static void reorderConnectedComponents(const Hypergraph& hg, ConnectedComponents& cc) {
		//compute new cc ids
		std::vector<nodeid> bucket_begin(hg.totalNodeWeight()+1, 0);
		for (auto csize : cc.ccsizes) {
			assert(csize > 0);
			assert(csize < hg.totalNodeWeight());	//don't call the DP with a connected hypergraph
			bucket_begin[csize]++;
		}
		std::partial_sum(bucket_begin.begin(), bucket_begin.end(), bucket_begin.begin());
		assert(bucket_begin.back() == cc.numComponents());
		std::vector<nodeid> old2new(cc.numComponents());
		for (nodeid c = 0; c < cc.numComponents(); c++) {
			nodeid new_ccid = bucket_begin[ cc.ccsizes[c]-1 ]++;
			assert(new_ccid < cc.numComponents());
			old2new[c] = new_ccid;
		}

		//apply new cc ids
		for (nodeid u = 0; u < hg.numNodes(); u++) { cc.node_component[u] = old2new[cc.node_component[u]]; }
		for (netid e = 0; e < hg.numHyperedges(); e++) { cc.hyperedge_component[e] = old2new[cc.hyperedge_component[e]]; }
		std::vector<nodeid> reordered_sizes(cc.numComponents());
		for (nodeid c_old = 0; c_old < cc.numComponents(); c_old++) { reordered_sizes[old2new[c_old]] = cc.ccsizes[c_old]; }
		std::swap(cc.ccsizes, reordered_sizes);
	}

	static void cutSingleHyperedgeConnectedComponent(nodeid nNodes, ParetoFront& pf) {
		for (nodeid u = 0; u < pf.cuts.size(); u++) {
			pf.cuts[u].smallerBlock = u;
			pf.cuts[u].largerBlock = nNodes - u;
			pf.cuts[u].nIsolatedNodes = 0;
			pf.cuts[u].cut = 1;
		}
		pf.cuts[0].cut = 0;
	}

	static void writeCutToFront(ParetoFront& result_front, std::vector<flow_t>& globalMinCut, nodeid numNodes, std::chrono::duration<double, std::micro> runtime) {
		nodeid global_perf_balance = numNodes / 2;
		for (nodeid u = 0; u <= global_perf_balance; u++) {
			if (globalMinCut[u] == MAX_FLOW) { continue; }
			CutInfo& resc = result_front.cuts[u];
			if (globalMinCut[u] < resc.cut) {
				resc.smallerBlock = u;
				resc.largerBlock = numNodes - u;
				resc.nIsolatedNodes = invalid_node;
				resc.cut = globalMinCut[u];
				resc.cutter_id = -2;
				resc.running_time_total = runtime;
			}
		}
	}

	//desired_balance_without_cut, first_cc_with_gap, ccsizes_at_first_gap
	static std::tuple<bool, nodeid, nodeid> find_gap(ConnectedComponents& cc, nodeid global_smb) {
		bool has_gap = false;
		bool desired_balance_without_cut = false;
		nodeid c = 0;
		nodeid ccsizes_sofar = 0;
		nodeid csize = 0;

		while (!has_gap && c < cc.numComponents() && !desired_balance_without_cut) {
			csize = cc.ccsizes[c];
			if (csize > ccsizes_sofar + 1) { has_gap = true; }
			else {
				if (ccsizes_sofar + csize >= global_smb) {
					desired_balance_without_cut = true;
				}
				ccsizes_sofar += csize;
				c++;
			}
		}
		if (!has_gap) {		//then we got perfect balance
			assert(desired_balance_without_cut);
		}
		return std::tie(desired_balance_without_cut, c, ccsizes_sofar);
	};

	static void includeNonGapCCs(nodeid global_perf_balance, nodeid numNodes, nodeid ccsizes_at_first_gap, std::vector<flow_t>& gmc, std::vector<flow_t>& ngmc) {
		ngmc = gmc;
		for (nodeid u = 1; u <= global_perf_balance; u++) {
			if (gmc[u] != MAX_FLOW) {
				nodeid reach_with_smaller_block = std::min(global_perf_balance, u + ccsizes_at_first_gap);
				ngmc[reach_with_smaller_block] = std::min(gmc[u], ngmc[reach_with_smaller_block]);

				nodeid n = numNodes - u;
				if (n <= global_perf_balance) {
					nodeid reach_with_larger_block = std::min(global_perf_balance, n + ccsizes_at_first_gap);
					ngmc[reach_with_larger_block] = std::min(gmc[u], ngmc[reach_with_larger_block]);
				}
			}
		}
		std::swap(ngmc, gmc);
	}

	/**
	 *
	 * @param cc the connected components
	 * @param ccsizes_at_first_gap  the number of nodes that are freely distributable
	 * @param first_cc_with_gap the first conneceted component with which not any balance is possible anymore
	 * @param numNodes
	 * @param global_result_front the result
	 * @param min_epsilon_min_blocksize = numNodes/2 = global_perf_balance if min eps = 0.0
	 * @return true, if min_epsilon_min_blocksize can be reached, false otherwise
	 */
	static bool subsetSum(ConnectedComponents& cc, nodeid ccsizes_at_first_gap, nodeid first_cc_with_gap, nodeid numNodes, std::vector<flow_t>& globalMinCut, nodeid min_epsilon_min_blocksize) {
		nodeid global_perf_balance = numNodes / 2;
		assert(globalMinCut.size() == global_perf_balance + 1);
		nodeid numRemainingNodes = numNodes - ccsizes_at_first_gap;
		std::vector<nodeid> index;
		std::vector<bool> summable(numRemainingNodes + 1, false);
		summable[0] = true;
		index.push_back(0);
		bool desired_balance_reached = false;

		for (nodeid c = first_cc_with_gap; c < cc.numComponents() && !desired_balance_reached; c++) {
			nodeid csize = cc.ccsizes[c];
			if (csize > global_perf_balance) { break; }	//this is the last one. csizes are sorted ascendingly
			std::vector<nodeid> insert_into_index;
			for (nodeid prev_summable : index) {
				nodeid e = prev_summable + csize;
				if (!summable[e]) {
					summable[e] = true;
					insert_into_index.push_back(e);
					//figure out if we can already get to the best of the requested balances. if so, stop early.
					//nodeid diff = std::min(e, ccsizes_at_first_gap);
					nodeid e2l = numRemainingNodes - e; //move e to larger block
					if ( 		(e <= min_epsilon_min_blocksize && e + ccsizes_at_first_gap >= min_epsilon_min_blocksize)
								||	(e2l <= min_epsilon_min_blocksize && e2l + ccsizes_at_first_gap >= min_epsilon_min_blocksize) )
					{
						desired_balance_reached = true;
						break;
					}
				}
			}
			index.insert(index.end(), insert_into_index.begin(), insert_into_index.end());
		}

		if (desired_balance_reached) {
			globalMinCut[min_epsilon_min_blocksize] = 0;
			return true;
		}

		nodeid max_achieved_balance = 0;
		auto write_range = [&](nodeid range_begin, nodeid inclusive_range_end) {
			if (range_begin > global_perf_balance) { return; }	//the for loop does not require it, but the success condition update does
			inclusive_range_end = std::min(global_perf_balance, inclusive_range_end);
			for (nodeid x = range_begin;  x <= inclusive_range_end; x++) {
				globalMinCut[x] = 0;
			}
			if (inclusive_range_end > max_achieved_balance) { max_achieved_balance = inclusive_range_end; }
		};
		std::sort(index.begin(), index.end());	//histogram sort would be possible, but why bother.
		for (std::size_t i = 1; i < index.size(); i++) {
			assert(summable[numRemainingNodes - index[i]] || numRemainingNodes - index[i] > global_perf_balance);
			assert(summable[numRemainingNodes - index[i-1]]  || numRemainingNodes - index[i-1] > global_perf_balance);
			assert(summable[index[i]]);
			assert(summable[index[i-1]]);
			//consider index[i-1] for smaller block. By symmetry, numRemainingNodes - index[i-1] can be summed to as well, and thus will appear later in the loop, which is equivalent to moving index[i-1] to the larger block
			//the range [index[i-1], index[i-1] + ccsizes_at_first_gap] can be achieved with zero cut.
			write_range(index[i-1], std::min(index[i-1] + ccsizes_at_first_gap, index[i] - 1));
		}
		write_range(index.back(), std::min(index.back() + ccsizes_at_first_gap, global_perf_balance));
		return max_achieved_balance >= min_epsilon_min_blocksize;
	};


	template<class Piercer, typename ExecuteAfterSTOption>
	static void run(const Hypergraph& global_hg, ParetoFront& result_front, ConnectedComponents& cc, State& state, std::vector<STOptions>& stOptions, ExecuteAfterSTOption easto)
	{
		reorderConnectedComponents(global_hg, cc);
		verifyInvariantsOfCCReordering(global_hg, cc, state);

		auto t = time_now();
		auto t_beginning = t;

		auto hgs_bound = ExtractConnectedComponents::extractConnectedComponentHypergraphs(global_hg, cc).first;
		const auto& hgs = hgs_bound;
		verifyInvariantsOfExtractedCCs(global_hg, cc, state, hgs);

		nodeid global_perf_balance = global_hg.totalNodeWeight() / 2;
		std::vector<flow_t> globalMinCutWithSubsetSum(global_perf_balance + 1, MAX_FLOW);
		globalMinCutWithSubsetSum[0] = 0;
		nodeid global_smb = state.getMaxSmallerBlocksize(global_hg.totalNodeWeight());
		assert(global_smb <= global_perf_balance);

		auto [no_gaps, first_cc_with_gap_bound, ccsizes_at_first_gap_bound] = find_gap(cc, global_smb);
		const auto & first_cc_with_gap = first_cc_with_gap_bound;
		const auto& ccsizes_at_first_gap = ccsizes_at_first_gap_bound;

		LOG << "No Gaps up to" << ccsizes_at_first_gap << V(first_cc_with_gap) << V(global_smb) << V(
				global_hg.totalNodeWeight());
		if (no_gaps) {
			assert(ccsizes_at_first_gap >= global_smb);
			for (nodeid u = 0; u <= std::min(ccsizes_at_first_gap, global_perf_balance); u++) { globalMinCutWithSubsetSum[u] = 0; }
			writeCutToFront(result_front, globalMinCutWithSubsetSum, global_hg.totalNodeWeight(), time_now() - t);
			STOptions stoption_finished_with_gaps;
			stoption_finished_with_gaps.description = "FinishedWithJustGaps";
			easto(stoption_finished_with_gaps, global_hg, result_front, state);
			return;
		}
		bool finish_with_subsetsum = subsetSum(cc, ccsizes_at_first_gap, first_cc_with_gap, global_hg.totalNodeWeight(), globalMinCutWithSubsetSum, global_smb);
		writeCutToFront(result_front, globalMinCutWithSubsetSum, global_hg.totalNodeWeight(), time_now() - t);
		if (finish_with_subsetsum) {
			STOptions stoption_finished_with_subsetsum;
			stoption_finished_with_subsetsum.description = "FinishedWithSubsetSum";
			easto(stoption_finished_with_subsetsum, global_hg, result_front, state);
			return;
		}

		bool split_only_one_component = first_cc_with_gap == cc.numComponents() - 1;
		if (split_only_one_component) LOG << "split only one component";
		assert(global_perf_balance > ccsizes_at_first_gap);	//otherwise find_gaps(..) would have been successful

		std::vector<ParetoFront> fronts;
		std::vector<EnsemblePartitioning> ensemblePartitionings;

		auto t_subsetsum = time_now();
		duration dur_subsetsum = duration(t_subsetsum - t);
		tlx::unused(dur_subsetsum);
		t = t_subsetsum;

		Hypergraph union_gap_cc_hgs;
		if ( std::any_of(stOptions.begin(), stOptions.end(), [](const STOptions& stoption){ return stoption.useFBT(); } ) ) {
			union_gap_cc_hgs = UnionDisjointHypergraphs::unionDisjointHypergraphs(hgs, [first_cc_with_gap=first_cc_with_gap](nodeid cc_id) { return cc_id >= first_cc_with_gap; });
			auto t_udhg = time_now();
			state.running_time_gen_fbt += duration(t_udhg - t); t = t_udhg;
		}

		flow_t external_cut_bound = MAX_FLOW;
		for (STOptions& stOption : stOptions) {
			//call patoh here to obtain partitions to be balanced
			Regrow::FBT fbp(global_hg.totalNodeWeight(), cc);
			fbp.dmco.first_cc_with_gap = first_cc_with_gap; fbp.dmco.cc_sizes_at_first_gap = ccsizes_at_first_gap;
			fbp.generateFBTWithPatoh(global_hg, union_gap_cc_hgs, result_front, stOption, state);

			external_cut_bound = std::min(external_cut_bound, result_front.getCut(global_smb).second.cut);
			for (auto& x : state.external_partitioner_results) {
				if (x.size_of_smaller_block >= global_smb) {
					external_cut_bound = std::min(external_cut_bound, x.cut);
				}
			}

			auto t_fbt = time_now(); state.running_time_gen_fbt += duration(t_fbt - t); t = t_fbt;

			auto cutComponents = [&]() {
				nodeid ccsizes_sofar = 0;
				std::vector<flow_t> gmc(global_perf_balance + 1, MAX_FLOW), ngmc(global_perf_balance + 1, MAX_FLOW);
				gmc[0] = 0;
				for (nodeid c = first_cc_with_gap; c < cc.numComponents(); c++) {
					LOG << "split component" << V(c);
					nodeid csize = cc.ccsizes[c];
					const Hypergraph& hg = hgs[c];
					LOG << "component hg n=" << hg.numNodes() << "m=" << hg.numHyperedges() << "p=" << hg.numPins()
						<< " nweight=" << hg.totalNodeWeight() << " mweight=" << hg.totalHyperedgeWeight() << " pweight=" << hg.totalNumPins();
					nodeid index_of_c = c - first_cc_with_gap;

					if (fronts.size() == index_of_c) {
						fronts.emplace_back(ParetoFront(hg.totalNodeWeight()));
						t = time_now();
						ensemblePartitionings.push_back(EnsemblePartitioning::buildPatohEnsemblePartitioning(hg, state));
						state.running_time_ensemble_partitioner += duration(time_now() - t);

					}

					t = time_now();

					ParetoFront& front = fronts[index_of_c];
					EnsemblePartitioning& ep = ensemblePartitionings[index_of_c];
					nodeid smb = hg.totalNodeWeight() / 2;
					if (split_only_one_component) {
						nodeid target_size_smaller_block = global_smb - ccsizes_at_first_gap;
						smb = std::min(hg.totalNodeWeight() / 2, target_size_smaller_block); //if we need fewer vertices than perfect balance, then only take as many as required
					}
					if (hg.numHyperedges() == 1) {
						cutSingleHyperedgeConnectedComponent(hg.totalNodeWeight(), front);
					}
					else {
						fbp.selectComponentForTerminalGeneration(c);
						MultiCutter::run<Piercer>(hg, front, ep, fbp, state, stOption, smb, external_cut_bound);
					}

					auto t_multicutter = time_now(); state.running_time_multi_cutter += duration(t_multicutter - t); t = t_multicutter;

					if (split_only_one_component) {
						for (nodeid u = 0; u < front.cuts.size(); u++) {
							if (front.cuts[u].cut != MAX_FLOW) {
								gmc[u] = std::min(gmc[u], front.cuts[u].cut);
								nodeid n = hg.totalNodeWeight() - u;
								if (n <= hg.totalNodeWeight()/2) {
									gmc[n] = std::min(gmc[u], front.cuts[u].cut);
								}
							}
						}
					}
					else {
						includeNextCCWithFlowCutterBasicVariantCacheFriendly(front, csize, ccsizes_sofar, global_perf_balance, ngmc, gmc);
					}
					ccsizes_sofar += csize;

					state.running_time_disconnected_dp += duration(time_now() - t);
				}

				t = time_now();
				includeNonGapCCs(global_perf_balance, global_hg.totalNodeWeight(), ccsizes_at_first_gap, gmc, ngmc);
				auto t_finished = time_now();
				state.running_time_disconnected_dp += duration(t_finished - t);

				duration total_runtime = t_finished - t_beginning;
				writeCutToFront(result_front, gmc, global_hg.totalNodeWeight(), total_runtime);
			};

			cutComponents();
			if (result_front.getCut(global_smb).second.cut == MAX_FLOW) {
				//no solution was possible (because terminals from FBT may contain more than n_c/2 nodes in block c)
				LOG << "enforce solvability";
				auto t_fbt_enforce = time_now();
				fbp.enforceSolvability(stOption, global_hg, union_gap_cc_hgs);
				state.running_time_gen_fbt += duration(time_now() - t_fbt_enforce);
				cutComponents();
			}
			easto(stOption, global_hg, result_front, state);
		}
	}

	static void includeNextCCWithFlowCutterBasicVariantCacheFriendlyWithoutIndex(ParetoFront& front, nodeid csize, nodeid ccsizes_sofar,
																				 nodeid global_perf_balance, std::vector<flow_t>& ngmc, std::vector<flow_t>& gmc) {
		nodeid perf_balance = csize / 2;
		for (nodeid i = 0; i <= std::min(ccsizes_sofar, global_perf_balance); i++) { ngmc[i] = gmc[i]; }
		for (nodeid i = 1; i <= std::min(ccsizes_sofar + csize, global_perf_balance); i++) {
			//move smaller block of C to smaller global block
			for (nodeid j = 0; j <= std::min(perf_balance, i); j++) {
				assert(i-j < gmc.size()); assert(j <= i); assert(j < front.cuts.size());
				if (front.cuts[j].isValid() && gmc[i-j] != MAX_FLOW) {
					if (ngmc[i] > gmc[i-j] + front.cuts[j].cut) { ngmc[i] = gmc[i-j] + front.cuts[j].cut; }
				}
			}
			//move larger block of C to smaller global block
			nodeid jmin = csize <= i ? 0 : csize - i;
			jmin = std::max<nodeid>(1, jmin);	//ignore j = 0
			for (nodeid j = jmin; j <= perf_balance; j++) {
				assert(i + j >= csize);
				if (front.cuts[j].isValid() && gmc[i + j - csize] != MAX_FLOW) {
					flow_t cut = gmc[i + j - csize] + front.cuts[j].cut;
					if (ngmc[i] > cut) { ngmc[i] = cut; }
				}
			}
		}
		std::swap(ngmc, gmc);
	}

	static void includeNextCCWithFlowCutterBasicVariantCacheFriendlyWithIndex(std::vector<nodeid>& valid_cuts, ParetoFront& front, nodeid csize, nodeid ccsizes_sofar,
																			  nodeid global_perf_balance, std::vector<flow_t>& ngmc, std::vector<flow_t>& gmc) {
		nodeid perf_balance = csize / 2;
		for (nodeid i = 0; i <= std::min(ccsizes_sofar, global_perf_balance); i++) { ngmc[i] = gmc[i]; }
		for (nodeid i = 1; i <= std::min(ccsizes_sofar + csize, global_perf_balance); i++) {
			//move smaller block of C to smaller global block
			nodeid maxJ = std::min(perf_balance, i);
			for (nodeid j : valid_cuts) {
				if (j > maxJ) { break; }
				assert(i - j < gmc.size());
				assert(j <= i);
				assert(j < front.cuts.size());
				if (gmc[i - j] != MAX_FLOW) {
					if (ngmc[i] > gmc[i - j] + front.cuts[j].cut) { ngmc[i] = gmc[i - j] + front.cuts[j].cut; }
				}
			}

			//move larger block of C to smaller global block
			nodeid jmin = csize <= i ? 0 : csize - i;
			jmin = std::max<nodeid>(1, jmin);	//ignore j = 0
			for (auto rit = valid_cuts.rbegin(); rit != valid_cuts.rend(); rit++) {
				nodeid j = *rit;
				if (j < jmin) { break; }
				assert(i + j >= csize);
				if (gmc[i + j - csize] != MAX_FLOW) {
					flow_t cut = gmc[i + j - csize] + front.cuts[j].cut;
					if (ngmc[i] > cut) { ngmc[i] = cut; }
				}
			}
		}
		std::swap(ngmc, gmc);
	}


	static void buildIndex(ParetoFront& front, std::vector<nodeid>& valid_cuts) {
		valid_cuts.clear();
		for (nodeid u = 0; u < front.cuts.size(); u++) {
			if (front.cuts[u].isValid()) {
				valid_cuts.push_back(u);
			}
		}
	}

	static bool useIndexForValidCuts(ParetoFront& front, std::vector<nodeid>& valid_cuts) {
		buildIndex(front, valid_cuts);
		if (valid_cuts.size() > front.cuts.size()/2) {
			valid_cuts.clear();
			return false;
		}
		return true;
	}

	static void includeNextCCWithFlowCutterBasicVariantCacheFriendly(ParetoFront& front, nodeid csize, nodeid ccsizes_sofar, nodeid global_perf_balance, std::vector<flow_t>& ngmc, std::vector<flow_t>& gmc) {
		std::vector<nodeid> valid_cuts;
		if (useIndexForValidCuts(front, valid_cuts)) {
			includeNextCCWithFlowCutterBasicVariantCacheFriendlyWithIndex(valid_cuts, front, csize, ccsizes_sofar, global_perf_balance, ngmc, gmc);
		}
		else {
			includeNextCCWithFlowCutterBasicVariantCacheFriendlyWithoutIndex(front, csize, ccsizes_sofar, global_perf_balance, ngmc, gmc);
		}

		{
			for (nodeid i = 0; i <= std::min(ccsizes_sofar + csize, global_perf_balance) && csize + i <= global_perf_balance; i++)
				if (ngmc[i] != MAX_FLOW)
					gmc[csize + i] = ngmc[i];

		}
	}

	static void verifyInvariantsOfCCReordering(const Hypergraph &global_hg, ConnectedComponents &cc, State &state) {
#ifndef NDEBUG
		assert(std::is_sorted(cc.ccsizes.begin(), cc.ccsizes.end()));
		std::vector<nodeid> ccsizes(global_hg.numNodes());
		for (nodeid u = 0; u < global_hg.numNodes(); u++) {
			ccsizes[cc.node_component[u]] += global_hg.nodeWeight(u);
		}
		for (nodeid c = 0; c < cc.numComponents(); c++) {
			assert(ccsizes[c] == cc.ccsizes[c]);
			assert(ccsizes[c] > 0);
		}
		for (nodeid c : cc.node_component) {
			assert(c <= global_hg.numNodes());
			assert(c < cc.numComponents());
		}
		for (nodeid c : cc.hyperedge_component) {
			assert(c <= global_hg.numNodes());
			assert(c < cc.numComponents());
		}
#endif
	}

	static void verifyInvariantsOfExtractedCCs(const Hypergraph& global_hg, ConnectedComponents& cc, State& state, const std::vector<Hypergraph>& hgs) {
#ifndef NDEBUG
		netid sumHyperedgeWeights = 0;
		nodeid sumNodeWeights = 0;
		nodeid totalNumPins = 0;
		for (nodeid c = 0; c < cc.numComponents(); c++) {
			assert(hgs[c].totalNodeWeight() == cc.ccsizes[c]);
			sumHyperedgeWeights += hgs[c].totalHyperedgeWeight();
			sumNodeWeights += hgs[c].totalNodeWeight();
			totalNumPins += hgs[c].totalNumPins();
		}
		assert(sumHyperedgeWeights == global_hg.totalHyperedgeWeight());
		assert(sumNodeWeights == global_hg.totalNodeWeight());
		assert(totalNumPins == global_hg.totalNumPins());
#endif
	}
};

}//namespace
