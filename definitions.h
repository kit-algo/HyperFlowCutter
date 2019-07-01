#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <random>
#include <algorithm>
#include "util/logger.h"
#include "util/timer.h"
#include <chrono>



template<typename int1, typename int2>
int1 ceil_div(int1 numerator, int2 denominator) {
	return numerator/denominator + (numerator % denominator == 0 ? 0 : 1);
}

template<typename id_t>
struct const_range {
	typename std::vector<id_t>::const_iterator __begin;
	typename std::vector<id_t>::const_iterator __end;
	inline typename std::vector<id_t>::const_iterator begin() const { return __begin; }
	inline typename std::vector<id_t>::const_iterator end() const { return __end; }
	const_range(typename std::vector<id_t>::const_iterator __begin, typename std::vector<id_t>::const_iterator __end) : __begin(__begin), __end(__end) { }
	const_range(std::vector<id_t>& vec, std::size_t a, std::size_t b) : const_range(vec.begin() + a, vec.begin() + b) {}
};


namespace aux {
	template<typename T> std::pair<T,T> minmax(T a, T b) {
		return a < b ? std::make_pair(a,b) : std::make_pair(b,a);
	}

	template<typename T>
	void min_to(T& a, T& b) { a = std::min(a,b); }
}

namespace hyper {
	using arcid = uint32_t;
	using nodeid = uint32_t;
	using label = nodeid;
	using netid = uint32_t;
	using partitionid = uint32_t;
	using flow_t = int32_t;

	//TODO substitute usages of nodeid and netid by Node and Hyperedge so the compiler catches more mistakes. also do NodeIndexable vectors
	struct Node {
		nodeid value;
	};

	struct Hyperedge {
		netid value;
	};

	static constexpr arcid invalid_arc = std::numeric_limits<arcid>::max();
	static constexpr nodeid invalid_node = std::numeric_limits<nodeid>::max();
	static constexpr netid invalid_hyperedge = std::numeric_limits<netid>::max();
	static constexpr label invalid_label = invalid_node;
	static constexpr partitionid invalid_block = invalid_node;

	static constexpr nodeid INF = std::numeric_limits<nodeid>::max() - 2;
	static constexpr flow_t MAX_FLOW = std::numeric_limits<flow_t>::max();

	struct ConnectedComponents {
		std::vector<nodeid> node_component;
		std::vector<nodeid> hyperedge_component;
		std::vector<nodeid> ccsizes;
		std::size_t numComponents() const { return ccsizes.size(); }

		bool isConnected() const {
			return numComponents() <= 1;
		}

		nodeid nodeComponent(nodeid u) const {
			return isConnected() ? 0 : node_component[u];
		}

		nodeid componentSize(nodeid c) const {
			return ccsizes[c];
		}

		static ConnectedComponents createCCOfConnectedHG(nodeid numNodes) {
			return { {}, {}, {numNodes} };
		}
	};

	struct ExternalPartitionerResult {
		flow_t cut;
		nodeid size_of_smaller_block;
	};

	enum MultiCutterExecution {
		Interleaved, Hybrid, Consecutive
	};

	struct PiercerOptions {
		bool use_threehop = true;
		bool use_random = true;
		bool avoid_augmenting_paths = true;
		bool use_ensemble = false;
		uint32_t numEnsemblePartitions = 0;
	};

	namespace Regrow {
		enum RegrowBlockStrategy {
			BisectStrat,
			EjectRandomFromBorderStrat
		};
		enum FBTPartitionStrategy {
			TwoWay,
			ThreeWay
		};
	}

	struct STOptions {
		uint32_t numRandomSt = 20;
		uint32_t numFarSt = 0;
		uint32_t numFurtherSt = 0;
		uint32_t numEnsembleSt = 0;

		uint32_t numDistinctFinishBalanceTerminals = 0;
		uint32_t numRepetitionsPerFinishBalanceTerminal = 1;
		uint32_t numExternalPartitionerCallsPerFinishBalanceTerminal = 1;
		double fbt_epsilon = 0.05;
		double fbt_fractional_blocksize = 0.5;
		std::string patoh_preset = "D";
		Regrow::RegrowBlockStrategy fbt_regrow_strat = Regrow::RegrowBlockStrategy::BisectStrat;
		Regrow::FBTPartitionStrategy fbt_partition_strat = Regrow::FBTPartitionStrategy::TwoWay;
		bool useFBT() const { return numDistinctFinishBalanceTerminals > 0 && numRepetitionsPerFinishBalanceTerminal > 0 && numExternalPartitionerCallsPerFinishBalanceTerminal > 0; }
		flow_t perfectly_balanced_cut_with_external_partitioner = MAX_FLOW;

		uint32_t getNumStPairs() { return numRandomSt + numFarSt + numFurtherSt + numEnsembleSt; }
		std::string description = "Default-ST-Options";
	};

	enum STMethod {
		Random, Far, Further, Ensemble, RandomDuplicatesAllowed, FBT
	};

	enum TypeOfFlowAlgorithm {
		Dinic, VertexDisjointEdmondsKarp
	};

	struct STPair {
		std::vector<nodeid> s,t;
		STMethod origin;
		std::size_t generatorID = 0;
		bool skip = false;
		bool s_is_coarse = false;
		bool t_is_coarse = false;
	};

	class Metrics {
	public:
		static nodeid largerBlockSize(const nodeid numNodes, const double eps) {
			nodeid perf_balance = numNodes/2 + (numNodes % 2);
			auto sm = static_cast<nodeid>( std::lround(std::floor((1.0+eps)/2.0 * numNodes)) );
			if (eps == 0.0) sm = perf_balance;
			return std::max(perf_balance, sm);
		}

		static nodeid smallerBlockSize(const nodeid numNodes, const double eps) {
			nodeid perf_balance = numNodes/2;
			auto sm = static_cast<nodeid>( std::lround(std::ceil((1.0-eps)/2.0 * numNodes)) );
			if (eps == 0.0) sm = perf_balance;
			return std::min(perf_balance, sm);
		}

		static double imbalance(const nodeid numNodes, const nodeid _smallerBlockSize) {
			double achieved_eps = 1.0 - (2.0 * static_cast<double>(_smallerBlockSize) / static_cast<double>(numNodes));
			if (_smallerBlockSize == numNodes / 2) { achieved_eps = 0.0; }
			return achieved_eps;
		}
	};

	class State {
	public:
		std::string hypergraphFile = "not loaded";
		std::string outputFilename = "cout";
		std::string stSpecificOutputFilename;

		std::string algo_name = "HyperFlowCutter";

		PiercerOptions piercerOptions;
		uint32_t seed = 855;

		MultiCutterExecution multiCutterExecution = MultiCutterExecution::Interleaved;

		TypeOfFlowAlgorithm type_of_flow_algorithm = TypeOfFlowAlgorithm::VertexDisjointEdmondsKarp;

		std::vector<double> epsilons = {0.2, 0.1, 0.03, 0.001, 0.0001, 0.00001, 0.0};

		int output_detail = 0;

		bool build_datastructures_during_grow_reachable = true;

		bool outputToCout = false;
		bool stSpecificOutput = false;
		bool onlyOutputTotalRunningTime = false;
		std::chrono::duration<double> totalRunningTime;

		duration running_time_ensemble_partitioner = duration(0.0);
		duration running_time_multi_cutter = duration(0.0);
		duration running_time_disconnected_dp = duration(0.0);
		duration running_time_gen_fbt = duration(0.0);
		std::vector<ExternalPartitionerResult> external_partitioner_results;

		uint32_t nDeletedHyperedges;
		nodeid maxHyperedgeSize;
		bool hypergraphIsConnectedAfterHyperedgeFiltering;

		nodeid getMaxSmallerBlocksize(nodeid n) {
			double min_eps = *std::min_element(epsilons.begin(), epsilons.end());
			return Metrics::smallerBlockSize(n, min_eps);
		}

		std::string getHypergraphName() { return hypergraphFile.substr(hypergraphFile.find_last_of('/') + 1); }

		uint32_t numSTPairs;
		uint32_t getNumStPairs() { return numSTPairs; }

	};
}

namespace LoggingInformation {
	void set_output_detail(int __output_detail);
	int get_output_detail();
}

namespace Random {
	uint32_t getSeed();
	void setSeed(uint32_t seed);
	std::mt19937& getRNG();

	template<typename T> T randomNumber(T a = 0, T b=std::numeric_limits<T>::max())
	{
		std::uniform_int_distribution<T> dist(a,b);
		return dist(getRNG());
	}

	template<typename T> std::vector<T> randomSequence(std::size_t n, T a = 0, T b=std::numeric_limits<T>::max())
	{
		std::vector<T> res;
		std::uniform_int_distribution<T> dist(a,b);
		for (std::size_t i = 0; i < n; i++) res.push_back(dist(getRNG()));
		return res;
	}

	template<typename T> const T& sampleOne(const std::vector<T>& seq) {
		return seq[ randomNumber<std::size_t>(0, seq.size()-1) ]; }

	template<typename T> std::pair<T,T> sampleTwoDisjoint(std::vector<T>& seq) {
		std::vector<T> sample_out(2);
		std::sample(seq.begin(), seq.end(), sample_out.begin(), 2, std::mt19937{randomNumber<std::size_t>()});
		return std::make_pair(sample_out[0], sample_out[1]);
	}

}
