#pragma once
#include <tlx/cmdline_parser.hpp>
#include "../definitions.h"
#include <iostream>
#include <stdexcept>

namespace hyper {

class ProcessCommandLineOptionsGeneralHyperFlowCutter {
public:
	static std::pair<STOptions, State> processCommandLineOptions(int argc, const char* argv[], std::string description = "HyperFlowCutter") {
		tlx::CmdlineParser cp;
		cp.set_description(description);

		State state;
		cp.add_param_string("hypergraph", state.hypergraphFile, "Hypergraph file");
		cp.add_param_string("output-table", state.outputFilename, "Where to write the output table. Set to cout to write to terminal");

		cp.add_uint('s', "seed", state.seed, "Random seed. Default 855");

		bool disable_aap = false, disable_random = false, disable_hop = false;
		cp.add_bool("disable-aap", disable_aap, "Disable avoiding augmenting paths during piercing.");
		cp.add_bool("disable-random", disable_random, "Disable randomized piercing decisions.");
		cp.add_bool("disable-hop", disable_hop, "Disable greedily minimizing mixed cut hyperedges.");
		bool disable_ensemble = false;
		cp.add_bool("disable-ensemble", disable_ensemble, "Disable ensemble piercing. Set if you want to use ensemble st-pairs without ensemble piercing.");
		cp.add_uint("ensemble", state.piercerOptions.numEnsemblePartitions, "# ensemble partitions used for piercing. Default 0. Currently only PaToH is supported as Ensemble Partitioner.");


		STOptions stOptions;
		cp.add_uint("st-random", stOptions.numRandomSt, "# random st pairs. Default 20.");
		cp.add_uint("st-further", stOptions.numFurtherSt, "# further pseudoperipheral st pairs. Default 0.");
		cp.add_uint("st-far", stOptions.numFarSt, "# far pseudoperipheral st pairs. Default 0.");
		cp.add_uint("st-ensemble", stOptions.numEnsembleSt, "# ensemble st pairs. Default 0. Requires --ensemble > 0.");

		std::string executor = "interleaved";
		cp.add_string("executor-type", executor, "determines how multiple st pair cutters are run. Options are interleaved, consecutive, hybrid. Default interleaved.");
		cp.add_bool("only-total-runningtime", state.onlyOutputTotalRunningTime, "Only output total running time. Default: off.");

		std::vector<std::string> string_epsilons;
		cp.add_stringlist("eps", string_epsilons, "Imbalances for which to write output. Default 0.2, 0.1, 0.03, 0.001, 0.0001, 0.00001, 0.0.");

		cp.add_string("vst-file", state.stSpecificOutputFilename, "Write output for every st-pair to this file. Set to cout to write to terminal.");

		cp.set_verbose_process(false);
		if (!cp.process(argc, argv)) { exit(-1); }
		//cp.print_result();

		state.numSTPairs = stOptions.getNumStPairs();

		state.piercerOptions.use_threehop = !disable_hop;
		state.piercerOptions.avoid_augmenting_paths = !disable_aap;
		state.piercerOptions.use_random = !disable_random;

		if (stOptions.numEnsembleSt > 0 && state.piercerOptions.numEnsemblePartitions == 0) { std::cerr << "--st-ensemble > 0 requires --ensemble > 0" << std::endl; exit(-1); }
		state.piercerOptions.use_ensemble = state.piercerOptions.numEnsemblePartitions > 0;

		if (disable_ensemble && state.piercerOptions.numEnsemblePartitions > 0) { state.piercerOptions.use_ensemble = false; }

		if (executor == "interleaved") { state.multiCutterExecution = MultiCutterExecution::Interleaved; }
		else if (executor == "consecutive") { state.multiCutterExecution = MultiCutterExecution::Consecutive; }
		else if (executor == "hybrid") { state.multiCutterExecution = MultiCutterExecution::Hybrid; }
		else { std::cerr << "Unknown executor-type option." << std::endl; exit(-1); }

		if (!string_epsilons.empty()) {
			state.epsilons.clear();
			for (auto& s : string_epsilons) {
				double x;
				try {
					x = std::stod(s);
				} catch(std::invalid_argument& ia) {
					std::cerr << "Imbalance parameter " << s << " could not be converted to double." << std::endl; exit(-1);
				}
				if (x < 0.0 || x > 1.0) {
					std::cerr << "Epsilon " << x << " out of range [0,1]." << std::endl; exit(-1);
				}
				state.epsilons.push_back(x);
			}
			std::sort(state.epsilons.begin(), state.epsilons.end(), std::greater<>());
		}

		if (state.outputFilename == "cout") { state.outputToCout = true; }
		if (!state.stSpecificOutputFilename.empty()) {
			state.stSpecificOutput = true;
			if (state.stSpecificOutputFilename == "cout") { state.stSpecificOutput = true; }
		}

		Random::setSeed(state.seed);

		return std::make_pair(stOptions,state);
	}
};

}//namespace