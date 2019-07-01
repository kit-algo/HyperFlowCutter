#pragma once

#include <string>
#include <fstream>
#include "../algorithm/multicutter.h"

namespace hyper {

	class CompatibleLineFormat {
	public:
		static void writeStringStream(State& state, std::stringstream& output) {
			output.sync();
			if (state.outputToCout) {
				std::cout << output.str() << std::flush;
			}
			else {
				std::ofstream f(state.outputFilename);
				f << output.str();
				f.close();
			}
		}

		static void writeSTOption(State& state, STOptions& sto, duration dur) {
			std::stringstream output;
			output << "HyperFlowCutter,";
			output << sto.description << "," << state.getHypergraphName() << "," << state.seed << "," << inSeconds(dur).count() << "\n";
			writeStringStream(state, output);
		}

		//Format: algoname,hg,eps,achieved_eps,seed,cut,runtime
		static void write(State& state, const ParetoFront& front) {
			std::stringstream output;
			std::string executor_desc;
			if (state.multiCutterExecution == MultiCutterExecution::Consecutive) { executor_desc = "Consecutive"; }
			if (state.multiCutterExecution == MultiCutterExecution::Interleaved) { executor_desc = "Interleaved"; }
			if (state.multiCutterExecution == MultiCutterExecution::Hybrid) { executor_desc = "Hybrid"; }
			//std::string algoname = "HyperFlowCutter-" + executor_desc + "-" + std::to_string(state.getNumStPairs());
			std::string algoname = "HyperFlowCutter";
			if (state.onlyOutputTotalRunningTime) {
				auto c = front.getCut(0.0).second;
				flow_t cut = c.cut + state.nDeletedHyperedges;
				output << algoname << "," << state.getHypergraphName() << "," << state.seed << "," << state.getNumStPairs() << "," << cut << "," << inSeconds(state.totalRunningTime).count() << "\n";
			}
			else {
				for (double eps : state.epsilons) {
					auto [smb, c] = front.getCut(eps);
					flow_t cut = c.cut + state.nDeletedHyperedges;
					double achieved_eps = Metrics::imbalance(front.nNodes, smb);
					output << algoname << "," << state.getHypergraphName() << "," << eps << "," << achieved_eps << "," << state.seed << "," << cut << "," << inSeconds(c.running_time_total).count() << "\n";
				}
			}
			writeStringStream(state, output);
		}

		static void writeEverySTPair(State& state, std::vector<ParetoFront>& fronts) {
			if (!state.stSpecificOutput) { return; }
			std::stringstream output;
			for (const ParetoFront& front : fronts) {
				for (double eps : state.epsilons) {
					auto b = front.getCut(eps);
					auto c = b.second;
					nodeid smb = b.first;
					flow_t cut = c.cut + state.nDeletedHyperedges;
					double achieved_eps = Metrics::imbalance(front.nNodes, smb);
					output << state.getHypergraphName() << "," << eps << "," << state.seed << "," << smb << "," << achieved_eps << "," << cut << "," << inSeconds(c.running_time_cutter).count() << "\n";
				}
			}
			std::swap(state.stSpecificOutputFilename, state.outputFilename);
			writeStringStream(state, output);
			std::swap(state.stSpecificOutputFilename, state.outputFilename);	//this hack sucks dude. just write clean code, already.
		}

		static std::string writeSTPair(State& state, const ParetoFront& front, STOptions& stOptions, uint16_t cutter_id) {
			std::stringstream output;
			std::stringstream algoname;
			algoname << "HyperFlowCutter_Executor-Consecutive_Piercing-";

			if (state.piercerOptions.use_threehop) { algoname << "3Hop"; } else { algoname << "No3Hop"; }
			if (state.piercerOptions.avoid_augmenting_paths) {algoname << "AAP";} else {algoname << "NoAAP";}
			if (state.piercerOptions.use_ensemble) {algoname << "Ensemble"; } else { algoname << "NoEnsemble"; }
			if (state.piercerOptions.use_random) {algoname << "Random";} else { algoname <<"NoRandom"; }
			algoname << "_ST-";

			if(stOptions.numRandomSt > 0) algoname << "Random";
			if (stOptions.numFurtherSt > 0) algoname << "Further";
			if (stOptions.numFarSt > 0) algoname << "Far";
			if (stOptions.numEnsembleSt > 0) algoname << "Ensemble";
			//auto epss = {0.0};
			for (double eps : state.epsilons) {
			//for ( double eps : epss) {
				auto b = front.getCut(eps);
				auto c = b.second;
				nodeid smb = b.first;

				flow_t cut = c.cut + state.nDeletedHyperedges;
				double achieved_eps = Metrics::imbalance(front.nNodes, smb);
				output << algoname.str() << "," << state.getHypergraphName() << "," << eps << "," << smb << "," << achieved_eps << "," << state.seed << "," << cutter_id << "," << cut << "," << inSeconds(c.running_time_cutter).count() << "\n";
			}
			return output.str();
		}
	};

	namespace HLA {
		/*
		class WriteSTPairs {
		public:
			static void write(const std::string& filename, const std::vector<STPair> &st) {
				std::ofstream f(filename);
				for (auto &x : st) {
					f << x.s << "\t" << x.t << "\t" << x.origin << "\n";
				}
				f << std::flush;
				f.close();
			}
		};
		*/

		struct Cut {
			bool is_placeholder;
			flow_t cut;
			double eps;
			double desired_eps;
			bool perf_bala;
			int cutter_id;
			duration runtime;
			duration runtime_of_stcutter;
		};

		class WriteCuts {
		public:
			static void write(const std::string& filename, const std::vector<Cut>& cuts) {
				std::ofstream f(filename);
				for (auto& x : cuts) {
					std::string bala = x.perf_bala ? "True" : "False";
					std::string placeholder = x.is_placeholder ? "True" : "False";
					auto dur = static_cast<second_duration>(x.runtime);
					auto dur_stcutter = static_cast<second_duration>(x.runtime_of_stcutter);
					f << placeholder << "\t" << x.cut << "\t" << x.eps << "\t" << x.desired_eps << "\t" << bala << "\t" << x.cutter_id << "\t" << dur.count() << "\t" << dur_stcutter.count() << "\n";
				}
				f << std::flush;
				f.close();
			}
		};

		class WriteSTCuts {
		public:
			static void write(const std::string& filename, const std::vector<std::vector<Cut>>& cuts) {
				std::ofstream f(filename);
				for (auto& y : cuts) {
					for (auto& x : y) {
						std::string bala = x.perf_bala ? "True" : "False";
						std::string placeholder = x.is_placeholder ? "True" : "False";
						auto dur_stcutter = static_cast<second_duration>(x.runtime_of_stcutter);
						f << placeholder << "\t" << x.cut << "\t" << x.eps << "\t" << x.desired_eps << "\t" << bala << "\t" << x.cutter_id << "\t" << dur_stcutter.count() << "\n";
					}
				}
				f << std::flush;
				f.close();
			}
		};


	}
}