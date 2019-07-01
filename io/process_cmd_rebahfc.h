#pragma once
#include <tlx/cmdline_parser.hpp>
#include "../definitions.h"
#include <iostream>
#include <stdexcept>

namespace hyper {
	class ProcessCMDHyPaToHFlowCutter {
	private:
		static Regrow::RegrowBlockStrategy parse_regrow_block_strategy(std::string& strat) {
			if (strat == "Bisect") return Regrow::RegrowBlockStrategy::BisectStrat;
			else if (strat == "EjectFromBorder") return Regrow::RegrowBlockStrategy::EjectRandomFromBorderStrat;
			else throw std::runtime_error("Unknown FBT regrow block strategy " + strat);
		}

		static Regrow::FBTPartitionStrategy parse_fbt_partition_strategy(std::string& strat) {
			if (strat == "2Way") return Regrow::FBTPartitionStrategy ::TwoWay;
			else if (strat == "3Way") return Regrow::FBTPartitionStrategy ::ThreeWay;
			else throw std::runtime_error("Unknown FBT partition strategy " + strat);
		}

		static TypeOfFlowAlgorithm parse_flow_algorithm(std::string& flow_algo) {
			if (flow_algo == "VertexDisjointEdmondsKarp" or flow_algo =="VDEK") return TypeOfFlowAlgorithm ::VertexDisjointEdmondsKarp;
			else if (flow_algo == "Dinic") return TypeOfFlowAlgorithm::Dinic;
			else throw std::runtime_error("Unknown flow algorithm " + flow_algo);
		}
	public:
		static std::pair<STOptions, State> processCommandLineOptions(int argc, const char* argv[]) {
			tlx::CmdlineParser cp;
			cp.set_description("Refinement_and_Balancing_HyperFlowCutter");

			State state;
			STOptions stOptions;
			stOptions.numRandomSt = 0;
			stOptions.description = "HyPaToHFlowCutter";

			state.hypergraphFile = "not loaded";
			cp.add_string('g', "hypergraph", state.hypergraphFile, "Hypergraph file");
			cp.add_uint('s', "seed", state.seed, "Random seed. Default 855");

			stOptions.numDistinctFinishBalanceTerminals = 1;
			cp.add_uint('d', "ndfbt", stOptions.numDistinctFinishBalanceTerminals, "Number of distinct FBTs generated. Default: 1");
			cp.add_uint('r', "nrep", stOptions.numRepetitionsPerFinishBalanceTerminal, "Number of HyperFlowCutter repetitions per distinct FBT. Default: 1");
			cp.add_uint('e', "nepcfbt", stOptions.numExternalPartitionerCallsPerFinishBalanceTerminal, "Number of external partitioner calls for generating one FBT. Default: 1");

			cp.add_double("fbt-eps", stOptions.fbt_epsilon, "Epsilon for external partitioner. Default 0.05");
			double eps = 0.0;
			cp.add_double("eps", eps, "Epsilon for HyPaToHFlowCutter. Default 0.");

			//cp.add_double("fbteps", stOptions.fbt_epsilon, "Balance parameter for FBT partitioning calls. Default 0.05");
			cp.add_double("fractionalblocksize", stOptions.fbt_fractional_blocksize, "Fractional max blocksize for FBT. Default 0.5");
			std::string str_regrow_strat = "Bisect";
			cp.add_string("regrow_strategy", str_regrow_strat, R"(FBT. How to shrink too large blocks. Options are "Bisect" (default) and "EjectFromBorder".)");
			std::string str_fbt_partition_strat = "2Way";
			cp.add_string("fbt_partition_strategy", str_fbt_partition_strat, R"(FBT. How to obtain terminals. Options are "2Way" (default) and "3Way")");

			std::string str_flow_algo = "Dinic";
			cp.add_string("flow-algo", str_flow_algo, R"(Flow algorithm. Options are "VertexDisjointEdmondsKarp" (default) (alias: "VDEK") and Dinic.)");

			cp.add_int("output-detail", state.output_detail, "Level of output detail. Default 0 (no output).");

			cp.add_string("patoh-preset", stOptions.patoh_preset, "PaToH preset. Q(uality), S(peed) or D(efault) (our default choice).");

			bool disable_build_datastructures_during_grow_reachable = false;
			cp.add_bool("disable-build-datastructures-during-grow-reachable", disable_build_datastructures_during_grow_reachable, "Build flow algorithm datastructures during growing reachable sides.");

			bool highly_detailed_output = false;
			cp.add_bool("--detailed-output", highly_detailed_output, "set to highest output detail");

			cp.set_verbose_process(false);
			if (!cp.process(argc, argv)) { exit(-1); }

			if (state.hypergraphFile == "not loaded") { throw std::runtime_error("No hypergraph given."); }

			state.build_datastructures_during_grow_reachable = !disable_build_datastructures_during_grow_reachable;

			state.type_of_flow_algorithm = parse_flow_algorithm(str_flow_algo);

			if (highly_detailed_output) state.output_detail = 10000;

			state.epsilons = { eps };

			stOptions.fbt_partition_strat = parse_fbt_partition_strategy(str_fbt_partition_strat);
			stOptions.fbt_regrow_strat = parse_regrow_block_strategy(str_regrow_strat);

			if (stOptions.fbt_fractional_blocksize > 0.5) { throw std::runtime_error("Fractional blocksize > 0.5"); }
			if (!stOptions.useFBT()) {
				throw std::runtime_error("FBT deactivated. Set all parameters --ndfbt/-d, --nrep/-r, --nepcfbt/-e to at least 1.");
			}

			std::stringstream algo_name;
			algo_name  	<< "ReBaHFC-" << stOptions.patoh_preset
						<< "-nrep=" << stOptions.numRepetitionsPerFinishBalanceTerminal
						<< "-ext_eps=" << stOptions.fbt_epsilon
						//<< "-eps=" << eps
						<< "-frac_size=" << stOptions.fbt_fractional_blocksize
						<< "-regrow=" << str_regrow_strat
					;
			state.algo_name = algo_name.str();


			return std::make_pair(stOptions, state);
		}
	};
}