#pragma once

#include "../definitions.h"
#include "finish_balance_partitioning.h"

namespace hyper {




class STChooser {

public:
	boost::dynamic_bitset<> is_terminal_node;
	std::vector<nodeid> potential_seed_nodes;
	std::size_t generatorID = 0;
	static constexpr bool debug = false;

	explicit STChooser(const Hypergraph &hg) : is_terminal_node(hg.numNodes()), potential_seed_nodes(
			hg.numNodes()) {
		std::iota(potential_seed_nodes.begin(), potential_seed_nodes.end(), 0);
	}

	void fillPotentialSeedNodes(const Hypergraph& hg) {
		potential_seed_nodes.clear();
		for (nodeid u = 0; u < hg.numNodes(); u++) {
			if (!is_terminal_node[u])
				potential_seed_nodes.push_back(u);
		}
	}

	std::vector<STPair> chooseRandomSTPairs(const Hypergraph &hg, uint64_t num_st_pairs) {
		generatorID++;
		fillPotentialSeedNodes(hg);
		std::vector<STPair> p;
		nodeid s, t;
		nodeid ind_s, ind_t;
		while (potential_seed_nodes.size() > 1 && p.size() < num_st_pairs) {
			std::uniform_int_distribution<nodeid> dist(0, static_cast<nodeid>(potential_seed_nodes.size()-1));
			do {
				ind_s = dist(Random::getRNG());
				ind_t = dist(Random::getRNG());
			} while (ind_s == ind_t);
			s = potential_seed_nodes[ind_s]; t = potential_seed_nodes[ind_t];
			is_terminal_node.set(s); is_terminal_node.set(t);
			potential_seed_nodes[ind_s] = potential_seed_nodes.back(); potential_seed_nodes.pop_back();
			potential_seed_nodes[ind_t] = potential_seed_nodes.back(); potential_seed_nodes.pop_back();
			p.push_back({ {s}, {t}, STMethod::Random, generatorID });
		}
		return p;
	}

	std::vector<STPair> chooseFarAwaySTPairs(const Hypergraph &hg, uint64_t num_st_pairs) {
		generatorID++;
		std::vector<STPair> p;
		nodeid s, t;
		LayeredQueue<nodeid> queue(hg.numNodes());
		boost::dynamic_bitset<> node_visited(hg.numNodes());
		boost::dynamic_bitset<> he_used(hg.numHyperedges());
		std::vector<nodeid> layer_ranges;
		for (uint64_t i = 0; i < num_st_pairs; i++) {
			t = invalid_node;
			while (!potential_seed_nodes.empty() && t == invalid_node) {
				do {
					std::uniform_int_distribution<nodeid> dist(0, static_cast<nodeid>(potential_seed_nodes.size()-1));
					nodeid ind = dist(Random::getRNG());
					s = potential_seed_nodes[ind];
					potential_seed_nodes[ind] = potential_seed_nodes.back();
					potential_seed_nodes.pop_back();
				} while (is_terminal_node[s] && !potential_seed_nodes.empty());
				if (is_terminal_node[s]) { return p; }	//couldn't select as many st-pairs as desired but at least return those found so far.
				LayeredBFS::run(hg, s, queue, node_visited, he_used, layer_ranges);
				if (layer_ranges.size() <= 2) { continue; }//it doesn't work with less than 2 layers. but less than three layers kinda sucks too. so try another seed node.
				nodeid layerend = layer_ranges.back(); layer_ranges.pop_back();
				nodeid layerfront = layer_ranges.back(); layer_ranges.pop_back();
				for (nodeid u : queue.range(layerfront, layerend)) {
					if (!is_terminal_node[u]) { t = u; break; }
				}
				if (t == invalid_node) {
					layerend = layerfront;
					layerfront = layer_ranges.back(); layer_ranges.pop_back();
					for (nodeid u: queue.range(layerfront, layerend)) {
						if (!is_terminal_node[u]) { t = u; break; }
					}
				}
			}
			if (t == invalid_node) { return p; }
			p.push_back({{s}, {t},STMethod::Far, generatorID});
			is_terminal_node.set(s); is_terminal_node.set(t);
		}
		return p;
	}

	std::vector<STPair> chooseFurtherAwaySTPair(const Hypergraph& hg, uint64_t num_st_pairs) {
		generatorID++;
		std::vector<STPair> p;
		nodeid s, t;
		LayeredQueue<nodeid> queue(hg.numNodes());
		boost::dynamic_bitset<> node_visited(hg.numNodes());
		boost::dynamic_bitset<> he_used(hg.numHyperedges());
		std::vector<nodeid> layer_ranges;
		while (!potential_seed_nodes.empty() && p.size() < num_st_pairs) {
			t = invalid_node;
			do {
				std::uniform_int_distribution<nodeid> dist(0, static_cast<nodeid>(potential_seed_nodes.size()-1));
				nodeid ind = dist(Random::getRNG());
				s = potential_seed_nodes[ind];
				potential_seed_nodes[ind] = potential_seed_nodes.back();
				potential_seed_nodes.pop_back();
			} while (is_terminal_node[s] && !potential_seed_nodes.empty());
			if (is_terminal_node[s]) { return p; }	//couldn't select as many st-pairs as desired but at least return those found so far.
			LayeredBFS::run(hg, s, queue, node_visited, he_used, layer_ranges);
			if (layer_ranges.size() <= 2) { continue; }//it doesn't work with less than 2 layers. but less than three layers kinda sucks too. so try another seed node.
			nodeid layerend = layer_ranges.back(); layer_ranges.pop_back();
			nodeid layerfront = layer_ranges.back(); layer_ranges.pop_back();
			for (nodeid u : queue.range(layerfront, layerend)) {
				if (!is_terminal_node[u]) { t = u; break; }
			}
			if (t == invalid_node) {
				layerend = layerfront;
				layerfront = layer_ranges.back(); layer_ranges.pop_back();
				for (nodeid u: queue.range(layerfront, layerend)) {
					if (!is_terminal_node[u]) { t = u; break; }
				}
			}
			if (t == invalid_node) continue;
			s = invalid_node;
			LayeredBFS::run(hg, t, queue, node_visited, he_used, layer_ranges);
			layerend = layer_ranges.back(); layer_ranges.pop_back(); layerfront = layer_ranges.back(); layer_ranges.pop_back();
			for (nodeid u : queue.range(layerfront, layerend)) {
				if (!is_terminal_node[u]) { s = u; break; }
			}
			if (s == invalid_node) {
				layerend = layerfront;
				layerfront = layer_ranges.back(); layer_ranges.pop_back();
				for (nodeid u: queue.range(layerfront, layerend)) {
					if (!is_terminal_node[u]) { s = u; break; }
				}
			}
			if (s == invalid_node) continue;
			is_terminal_node.set(s); is_terminal_node.set(t);
			p.push_back({{s},{t},STMethod::Further, generatorID});
		}
		return p;
	}

	std::vector<STPair> chooseEnsembleSTPairs(const Hypergraph& hg, State& state, STOptions& stOptions, EnsemblePartitioning& ep) {
		std::vector<STPair> res;
		if (stOptions.numEnsembleSt == 0) { return res; }

		auto subsets = ep.intersectionPartition.membersOfAllSubsets().second;
		std::sort(subsets.begin(), subsets.end(), [](const auto& lhs, const auto& rhs) { return lhs.size() < rhs.size(); });
		while (!subsets.empty() && subsets.back().size() > hg.numNodes() / 2) {
			subsets.back().pop_back();
		}

		if (subsets.size() < 2) {
			std::cout << state.getHypergraphName() << "," << state.seed << ", too few blocks in ensemble intersection partition." << std::endl;
			exit(-1);
		}

		for (std::size_t i = 0; i < ep.numUsedEnsembleSTPairs && subsets.size() >= 2; i++) {
			subsets.pop_back();
			subsets.pop_back();
		}

		//always match the two largest
		while (res.size() < stOptions.numEnsembleSt && subsets.size() >= 2) {
			//give every ensemble pair a new generator ID to indicate that it is not possible to use the same contracted hypergraph for them.
			res.push_back({subsets[subsets.size()-2], subsets[subsets.size()-1], STMethod::Ensemble, ++generatorID});
			subsets.pop_back();
			subsets.pop_back();
			ep.numUsedEnsembleSTPairs++;
		}
		return res;
	}

	std::vector<STPair> chooseFinishBalanceTerminals(const Hypergraph& hg, State& state, STOptions& stOptions, Regrow::FBT& fbp) {
		return fbp.generateTerminalPairs(stOptions, ++generatorID);
	}


	void fillWithDuplicatesAllowed(const Hypergraph& hg, uint64_t num, std::vector<STPair>& pairs) {
		nodeid s = invalid_node, t = invalid_node;
		std::uniform_int_distribution<nodeid> dist(0, static_cast<nodeid>(hg.numNodes()-1));
		for (uint64_t i = 0; i < num; i++) {
			do {
				if (hg.numNodes() == 2) { s = 0; t = 1; }
				else {
					s = dist(Random::getRNG()); t = dist(Random::getRNG());
				}
			} while (s ==t);
			pairs.push_back({{s}, {t},STMethod::RandomDuplicatesAllowed, generatorID});	//generatorID not incremented. they belong to the random terminals!
		}
	}

	std::tuple<std::vector<STPair>, uint64_t, uint64_t, uint64_t> getSTPairs(const Hypergraph& hg, State& state,
																			 STOptions& stOptions, EnsemblePartitioning& ep, Regrow::FBT& fbp) {
		uint64_t numFurtherAwaySTPairs = stOptions.numFurtherSt, numFarAwaySTPairs= stOptions.numFarSt, numRandomSTPairs = stOptions.numRandomSt;
		std::vector<STPair> p0 = chooseEnsembleSTPairs(hg, state, stOptions, ep);
		numRandomSTPairs += (stOptions.numEnsembleSt - p0.size());
		//LOG << "further";
		std::vector<STPair> p1 = chooseFurtherAwaySTPair(hg, numFurtherAwaySTPairs);
		fillPotentialSeedNodes(hg);
		numFarAwaySTPairs += (numFurtherAwaySTPairs - p1.size());
		numFurtherAwaySTPairs = p1.size();
		//LOG << V(numFurtherAwaySTPairs) << "far";
		std::vector<STPair> p2 = chooseFarAwaySTPairs(hg, numFarAwaySTPairs);
		//LOG << V(numFarAwaySTPairs) << " random";
		numRandomSTPairs += (numFarAwaySTPairs - p2.size());
		numFarAwaySTPairs = p2.size();
		fillPotentialSeedNodes(hg);
		std::vector<STPair> p3 = chooseRandomSTPairs(hg, numRandomSTPairs);
		if (p3.size() != numRandomSTPairs) {
			fillWithDuplicatesAllowed(hg, numRandomSTPairs - p3.size(), p3);
		}
		std::vector<STPair> p4 = chooseFinishBalanceTerminals(hg, state, stOptions, fbp);
		p0.insert(p0.end(), p1.begin(), p1.end());
		p0.insert(p0.end(), p2.begin(), p2.end());
		p0.insert(p0.end(), p3.begin(), p3.end());
		p0.insert(p0.end(), p4.begin(), p4.end());
		return std::tie(p0, numFurtherAwaySTPairs, numFarAwaySTPairs, numRandomSTPairs);
	}
};

}//namespace