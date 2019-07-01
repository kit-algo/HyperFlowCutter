#pragma once

#include <tlx/cmdline_parser.hpp>
#include "../definitions.h"
#include <iostream>
#include <stdexcept>

namespace hyper {
	class ProcessCMDStaggeredHyperFlowCutter{
	private:
		static TypeOfFlowAlgorithm parse_flow_algorithm(std::string& flow_algo) {
			if (flow_algo == "VertexDisjointEdmondsKarp" or flow_algo =="VDEK") return TypeOfFlowAlgorithm ::VertexDisjointEdmondsKarp;
			else if (flow_algo == "Dinic") return TypeOfFlowAlgorithm::Dinic;
			else throw std::runtime_error("Unknown flow algorithm " + flow_algo);
		}
	public:
		static State processCommandLineOptions(int argc, const char* argv[]) {
			tlx::CmdlineParser cp;
			cp.set_description("StaggeredHyperFlowCutter");

			State state;
			state.hypergraphFile = "not loaded";
			cp.add_param_string("hypergraph", state.hypergraphFile, "Hypergraph file");
			cp.add_param_uint("seed", state.seed, "Random seed. Default 855");
			std::string str_flow_algo = "VertexDisjointEdmondsKarp";
			cp.add_string("flow-algo", str_flow_algo, R"(Flow algorithm. Options are "VertexDisjointEdmondsKarp" (default) (alias: "VDEK") and Dinic.)");
			cp.add_int("output-detail", state.output_detail, "Level of output detail. Default 0 (no output).");

			bool disable_build_datastructures_during_grow_reachable = false;
			cp.add_bool("disable-build-datastructures-during-grow-reachable", disable_build_datastructures_during_grow_reachable, "Build flow algorithm datastructures during growing reachable sides.");

			bool highly_detailed_output = false;
			cp.add_bool("--detailed-output", highly_detailed_output, "set to highest output detail");

			cp.set_verbose_process(false);
			if (!cp.process(argc, argv)) { exit(-1); }
			if (state.hypergraphFile == "not loaded") { throw std::runtime_error("No hypergraph given."); }
			state.type_of_flow_algorithm = parse_flow_algorithm(str_flow_algo);

			state.build_datastructures_during_grow_reachable = !disable_build_datastructures_during_grow_reachable;

			if (highly_detailed_output) state.output_detail = 10000;
			state.epsilons = {0.0};
			state.piercerOptions.numEnsemblePartitions = 10;
			//remaining piercer options get ignored! rightfully so
			state.outputToCout = true;
			state.multiCutterExecution = MultiCutterExecution::Interleaved;
			return state;
		}
	};
}//namespace