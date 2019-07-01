#pragma once

#include "../datastructure/hypergraph.h"
#include "../datastructure/partition.h"
#include "patoh.h"



using hyper::Partition;
using hyper::Hypergraph;
using hyper::nodeid;

class PaToHInterface {
public:

	struct OurPatohParameters {
		explicit OurPatohParameters(double epsilon=0.0) : use_target_weights(false), target_weights(), epsilon(epsilon), k(2) {}
		bool use_target_weights;
		std::vector<float> target_weights;
		double epsilon;
		int k;
	};

	static std::vector<Partition> threeWayPartitionWithTwoEqualSizedBlocksWithPatoh(const Hypergraph& hg, std::vector<int>& seeds,
																					double fractional_size_of_equal_sized_blocks=0.45,
																					double epsilon=0.05,
																					std::string preset = "D") {
		if (fractional_size_of_equal_sized_blocks > 0.5) {
			throw std::runtime_error("3-way partition with two equal sized blocks which should both have 0.5*n is not possible.");
		}
		OurPatohParameters p;
		p.use_target_weights = true;
		p.k = 3;
		p.epsilon = epsilon;
		p.target_weights = { float(fractional_size_of_equal_sized_blocks), float(fractional_size_of_equal_sized_blocks), float(1.0 - 2.0*fractional_size_of_equal_sized_blocks)};
		return runPatoh(hg, seeds, p, preset);
	}

	static std::vector<Partition> bipartitionWithMaxSmallerBlocksizeWithPatoh(const Hypergraph& hg,
																			  std::vector<int>& seeds,
																			  nodeid max_smaller_blocksize,
																			  double epsilon=0.05,
																			  std::string preset = "D") {
		OurPatohParameters p;
		p.use_target_weights = true;
		p.k = 2;
		p.epsilon = epsilon;
		p.target_weights = { float(max_smaller_blocksize), float(hg.totalNodeWeight() - max_smaller_blocksize) };
		return runPatoh(hg, seeds, p, preset);
	}


	static std::vector<Partition> bisectWithPatoh(const Hypergraph &hg,
												  std::vector<int> &seeds,
												  double epsilon=0.0,
												  std::string preset = "D") {
		return runPatoh(hg, seeds, OurPatohParameters(epsilon), preset);
	}


	static std::vector<Partition> runPatoh(const Hypergraph& hg, std::vector<int>& seeds, OurPatohParameters params, std::string str_preset = "D") {
		if (hg.hasCoarseNodes()) throw std::runtime_error("At the moment we do not support external partitioners with coarse hypergraphs. It is implemented, but not tested.");
		if (params.use_target_weights && params.target_weights.size() != static_cast<std::size_t>(params.k)) {
			throw std::runtime_error("Number of specified targets weights do not match k. k=" + std::to_string(params.k) + " . target_weights.size()=" + std::to_string(params.target_weights.size()));
		}
		std::vector<Partition> result(seeds.size(), Partition(hg.numNodes()));
		if (seeds.empty()) { return result; }

		std::vector<int> vec_partition(hg.numNodes());
		std::vector<int> vec_partweights(params.k, 0);
		std::vector<int> vec_cwghts(hg.numNodes(), 1);

		std::vector<int> vec_pins(hg.numPins());
		std::vector<int> vec_xpins(hg.numHyperedges() + 1);
		{
			//copy operator does not work
			auto &hgpins = hg.getPins();
			for (std::size_t i = 0; i < hg.numPins(); i++) { vec_pins[i] = static_cast<int>(hgpins[i]); }
			auto &hefo = hg.getHyperedgeFirstOut();
			for (std::size_t i = 0; i < hefo.size(); i++) { vec_xpins[i] = static_cast<int>(hefo[i]); }
			for (nodeid u = 0; u < hg.numNodes(); u++) { vec_cwghts[u] = hg.nodeWeight(u); }
		}

		PaToH_Parameters args;
		int preset = -1;
		if (str_preset == "D") { preset = PATOH_SUGPARAM_DEFAULT; }
		else if (str_preset == "Q") { preset = PATOH_SUGPARAM_QUALITY; }
		else if (str_preset == "S") { preset = PATOH_SUGPARAM_SPEED; }
		else { throw std::runtime_error("Unknown PaToH preset" + str_preset); }

		PaToH_Initialize_Parameters(&args, PATOH_CUTPART, preset);
		args._k = params.k;
		args.doinitperm = 0;
		args.final_imbal = params.epsilon;

		if (str_preset == "Q") {
			args.MemMul_Pins = 100;
			args.MemMul_CellNet = 100;
			args.MemMul_General = 100;
		}

		/*
		std::cout << "coarsener " << args.crs_alg << " refiner " << args.ref_alg << std::endl;

		std::cout << "init imbal " << args.init_imbal << " final imbal " << args.final_imbal << std::endl;

		std::cout << "fixed netsize thrsh " << args.bisec_fixednetsizetrsh << " netsize thrsh " << args.bisec_netsizetrsh << " netsize part mul thrsh " << args.bisec_partmultnetsizetrsh << std::endl;

		std::cout << "visit order " << args.crs_VisitOrder << std::endl;

		std::cout << "bigV " << args.bigVcycle << " smallV " << args.smallVcycle << std::endl;

		std::cout << "ref max neg move " << args.ref_maxnegmove << " dyn lock " << args.ref_dynamiclockcnt << std::endl;


		*/
		args.outputdetail = PATOH_OD_HIGH;
		//args.outputdetail = PATOH_OD_LOW;

		int c, n, nconst, *cwghts, *xpins, *pins, *partvec, cut, *partweights;
		float* targetweights;
		n = hg.numHyperedges();
		c = hg.numNodes();
		nconst = 1;
		partvec = vec_partition.data();
		partweights = vec_partweights.data();
		cwghts = vec_cwghts.data();
		pins = vec_pins.data();
		xpins = vec_xpins.data();
		if (params.use_target_weights) {
			targetweights = params.target_weights.data();
		}
		else {
			targetweights = NULL;
		}

		PaToH_Alloc(&args, c, n, nconst, cwghts, NULL, xpins, pins);
		for (std::size_t i = 0; i < seeds.size(); i++) {
			args.seed = seeds[i];
			PaToH_Part(&args, c, n, nconst, 0, cwghts, NULL, xpins, pins, targetweights, partvec, partweights, &cut);
			for (hyper::nodeid u = 0; u < hg.numNodes(); u++) { result[i][u] = static_cast<hyper::nodeid>(vec_partition[u]); }
		}

		PaToH_Free();
		return result;
	}

};
