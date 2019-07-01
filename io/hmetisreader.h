#pragma once

#include <vector>
#include "../definitions.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "../datastructure/hypergraph.h"

namespace hyper {
	enum class HYPERGRAPH_TYPE : uint8_t {
		UNWEIGHTED = 0,
		EDGEWEIGHTS = 1,
		NODEWEIGHTS = 10,
		EDGEANDNODEWEIGHTS = 11,
	};

	struct hypergraph_info {
		HYPERGRAPH_TYPE type;
		nodeid num_nodes;
		netid num_nets;
		std::vector<nodeid> netsizes;
		std::vector<nodeid> pins;
		netid num_oversized_nets;
	};

	inline void mgetline(std::ifstream& f, std::string& line) {
		std::getline(f, line);
		while (line[0] == '%') {
			std::getline(f,line);
		}
	}

	void read_hmetis(std::string& filename, hypergraph_info& hg);

	void read_hmetis_with_hyperedge_filtering(std::string& filename, hypergraph_info& hg);

	class HMetisReader {
	public:
		static Hypergraph read(std::string& filename) {
			hypergraph_info hgi;
			read_hmetis(filename, hgi);
			return Hypergraph(hgi.num_nodes, hgi.netsizes, std::move(hgi.pins));
		}

		static Hypergraph readWithHyperedgeFiltering(State& state) {
			hypergraph_info hgi;
			read_hmetis(state.hypergraphFile, hgi);
			double max_eps = *std::max_element(state.epsilons.begin(), state.epsilons.end());
			state.maxHyperedgeSize = Metrics::largerBlockSize(hgi.num_nodes, max_eps);
			return Hypergraph::buildWithoutLargeHyperedges(hgi.num_nodes, hgi.netsizes, std::move(hgi.pins), state.maxHyperedgeSize, state.nDeletedHyperedges);
		}
	};
}//namespace hyper
