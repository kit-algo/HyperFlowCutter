#include "io/process_cmd_staggered_hyperflowcutter.h"
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

		auto easto = [&](STOptions& sto, const Hypergraph& hg, ParetoFront& front, State& state) {
			CompatibleLineFormat::writeSTOption(state, sto, time_now() - t);
			CompatibleLineFormat::write(state, front);
			//distinguishable runtime format. all in seconds: ensemble, dp, gen_fbt, multicutter
			std::cout
					<< "runtime ensemble: " << inSeconds(state.running_time_ensemble_partitioner).count() << "s "
					<< "runtime dp: " << inSeconds(state.running_time_disconnected_dp).count() << "s "
					<< "runtime gen fbt: " << inSeconds(state.running_time_gen_fbt).count() << "s "
	 				<< "runtime multicutter: " << inSeconds(state.running_time_multi_cutter).count() << "s"
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

	std::vector<STOptions> shfc1_test_config(State& state) {
		STOptions sto;
		sto.numRandomSt = 1;
		sto.numEnsembleSt = 0;
		sto.description = "FirstWave1Random";

		std::vector<STOptions> res;
		res.push_back(sto);

		state.piercerOptions.numEnsemblePartitions = 10;
		//remaining piercer options get ignored! rightfully so
		state.outputToCout = true;
		state.numSTPairs = 1;
		state.multiCutterExecution = MultiCutterExecution::Interleaved;
		return res;
	}


	std::vector<STOptions> shfc100_expconfig_stoptions(State& state) {
		STOptions bigchunk;
		bigchunk.numRandomSt = 78;
		bigchunk.numEnsembleSt = 2;
		bigchunk.description = "FourthWave78Random2Ensemble";

		STOptions thirdgroup;
		thirdgroup.numRandomSt = 14;
		thirdgroup.description = "ThirdWave14Random";

		STOptions smallgroup;
		smallgroup.numRandomSt = 5;
		smallgroup.description = "SecondWave5Random";

		STOptions ensemble;
		ensemble.numRandomSt = 0;
		ensemble.numEnsembleSt = 1;
		ensemble.description = "FirstWave1Ensemble";

		std::vector<STOptions> res;
		res.push_back(ensemble);
		res.push_back(smallgroup);
		res.push_back(thirdgroup);
		res.push_back(bigchunk);

		state.numSTPairs = 100;
		return res;
	}

	void disable_ensemble_if_not_used(State& state, std::vector<STOptions>& stOptions) {
		//in StaggeredHyperFlowCutter we don't use ensemble piercing, since it did improve solution quality.
		if (std::all_of(stOptions.begin(), stOptions.end(), [](const auto& sto) { return sto.numEnsembleSt == 0; })) {
			state.piercerOptions.numEnsemblePartitions = 0;
		}
	}
}


int main(int argc, const char* argv[]) {
	auto state = hyper::ProcessCMDStaggeredHyperFlowCutter::processCommandLineOptions(argc, argv);
	//auto stOptions = hyper::shfc1_test_config(state);
	auto stOptions = hyper::shfc100_expconfig_stoptions(state);
	hyper::disable_ensemble_if_not_used(state, stOptions);
	hyper::hidden_main(state, stOptions);

}
