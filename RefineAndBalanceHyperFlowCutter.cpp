#include "io/process_cmd_rebahfc.h"
#include "algorithm/multicutter_disconnected.h"
#include "io/hlafc_output.h"
#include "algorithm/hla_piercing_heuristics/hla_basic_piercing_filters.h"

namespace hyper {


	void hidden_main(State& state, std::vector<STOptions>& stOptions) {
		Random::setSeed(state.seed);
		LoggingInformation::set_output_detail(state.output_detail);
		Hypergraph hg = HMetisReader::readWithHyperedgeFiltering(state);
		ConnectedComponents cc = ConnectivityCheck::connectedComponents(hg);
		state.hypergraphIsConnectedAfterHyperedgeFiltering = cc.numComponents() <= 1;
		ParetoFront front(hg.totalNodeWeight());
		auto t = time_now();
		auto easto = [&](STOptions& sto, const Hypergraph& hg, ParetoFront& lambda_front, State& state) {
			//cuts of external partitioner and runtime can be accessed via state.external_partitioner_results and state.running_time_gen_fbt

			double desired_eps = state.epsilons[0];
			auto desired_smb = Metrics::smallerBlockSize(hg.totalNodeWeight(), desired_eps);

			auto [achieved_smb, c] = lambda_front.getCut(desired_eps);
			flow_t hfc_cut = c.cut + state.nDeletedHyperedges;
			double achieved_eps = Metrics::imbalance(hg.totalNodeWeight(), achieved_smb);

			//Format: algoname,hg,eps,achieved_eps,seed,cut,runtime
			std::cout
				<< state.algo_name << "," << state.getHypergraphName() << ","
				<< desired_eps << "," << achieved_eps << ","
				<< state.seed << "," << hfc_cut << ","
				<< inSeconds(duration(time_now() - t)).count() << std::endl;


			flow_t external_cut = MAX_FLOW;
			for (auto& x : state.external_partitioner_results) {
				if (x.size_of_smaller_block >= desired_smb) {
					external_cut = std::min(external_cut, x.cut);
				}
			}

			std::cout << "ExternalCut," << external_cut << std::endl;

			//distinguishable runtime format. all in seconds: ensemble, dp, gen_fbt, multicutter
			//std::cout << "runtime ensemble" << "runtime dp" << "runtime gen fbt" << "runtime multicutter" << std::endl;
			std::cout
				<< inSeconds(state.running_time_ensemble_partitioner).count() << ","
				<< inSeconds(state.running_time_disconnected_dp).count() << ","
				<< inSeconds(state.running_time_gen_fbt).count() << ","
				<< inSeconds(state.running_time_multi_cutter).count()
				<< std::endl;
		};

		if (cc.numComponents() <= 1) {
			auto t_ensemble = time_now();
			EnsemblePartitioning ep = EnsemblePartitioning::buildPatohEnsemblePartitioning(hg, state);
			state.running_time_ensemble_partitioner += duration(time_now() - t_ensemble);
			nodeid smb = state.getMaxSmallerBlocksize(hg.totalNodeWeight());
			MultiCutter::run<IgnoreOptionsPiercing::PiercingPipeline>(hg, front, state, ep, stOptions, smb, easto);
		}
		else {
			DisconnectedMultiCutter::run<IgnoreOptionsPiercing::PiercingPipeline>(hg, front, cc, state, stOptions, easto);
		}
	}
}

int main(int argc, const char* argv[]) {
	auto [stOption, state] = hyper::ProcessCMDHyPaToHFlowCutter::processCommandLineOptions(argc, argv);
	std::vector<hyper::STOptions> stOptions(stOption.numDistinctFinishBalanceTerminals, stOption);
	hyper::hidden_main(state, stOptions);
}
