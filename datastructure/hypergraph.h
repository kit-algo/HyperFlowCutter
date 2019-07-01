#pragma once

#include <cassert>
#include <numeric>
#include "../definitions.h"
#include <vector>
#include <random>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

namespace hyper {
	class HypergraphWithNodeWeightsWithoutHyperedgeWeights {
	private:
		nodeid num_contracted_nodes = 0;			// = sum( coarse_node_weights ) - numCoarseNodes()
		netid num_single_pin_hyperedges = 0;
		std::size_t num_contracted_pins = 0;
		std::vector<nodeid> coarse_node_weights;

		std::vector<netid> node_first_out;
		std::vector<netid> incident_hyperedges;
		std::vector<nodeid> hyperedge_first_out;
		std::vector<nodeid> pins;

		void init(std::vector<nodeid>& hyperedge_sizes) {
			if (totalNodeWeight() == 0) throw std::runtime_error("Cannot create hypergraph with zero total vertex weight.");
			if (numNodes() == 0) throw std::runtime_error("Cannot create hypergraph with zero vertices.");

			for (nodeid pin : pins) { node_first_out[pin+1]++; }
			std::partial_sum(node_first_out.begin(), node_first_out.end(), node_first_out.begin());
			for (netid e = 1; e < numHyperedges() + 1; e++) {
				hyperedge_first_out[e] = hyperedge_first_out[e-1] + hyperedge_sizes[e-1];
				for (nodeid pinind = hyperedge_first_out[e-1]; pinind < hyperedge_first_out[e]; pinind++) {
					nodeid pin = pins[pinind];
					incident_hyperedges[node_first_out[pin]++] = e-1;
				}
			}
			for (nodeid u = numNodes()-1; u > 0; u--) { node_first_out[u] = node_first_out[u-1]; }
			node_first_out[0] = 0;
		}

	public:
		using const_pin_iterator = std::vector<nodeid>::const_iterator;
		using pin_iterator = const_pin_iterator;
		using pin_range = const_range<nodeid>;
		using const_hyperedge_iterator = std::vector<netid>::const_iterator;
		using hyperedge_iterator = const_hyperedge_iterator;
		using hyperedge_range = const_range<netid>;


		/**
		 * empty hypergraph
		 */
		HypergraphWithNodeWeightsWithoutHyperedgeWeights() = default;


		/**
		 * Standard hypergraph without any contractions
		 * @param num_nodes	 = numIndexableNodes
		 * @param hyperedge_sizes
		 * @param _pins
		 */
		HypergraphWithNodeWeightsWithoutHyperedgeWeights(nodeid num_nodes, std::vector<nodeid>& hyperedge_sizes, std::vector<nodeid> _pins) :
				node_first_out(num_nodes+1, 0),
			 	incident_hyperedges(_pins.size()),
				hyperedge_first_out(hyperedge_sizes.size()+1, 0),
				pins(std::move(_pins))
		{
			init(hyperedge_sizes);
		}

		/**
		 * Hypergraph with every super_node contracted into one node. Single pin hyperedges are eliminated. Parallel hyperedges are not eliminated, since we do not support weighted hyperedges.
		 * @param fine_hg
		 * @param coarse_nodes for each coarse node, the collection of its fine nodes in fine_hg
		 */
		HypergraphWithNodeWeightsWithoutHyperedgeWeights(
				const HypergraphWithNodeWeightsWithoutHyperedgeWeights& fine_hg,
				const std::vector<std::vector<nodeid>>& coarse_nodes):
						num_contracted_nodes(
								std::accumulate(coarse_nodes.cbegin(), coarse_nodes.cend(), nodeid(0), [](const auto& n_contr_nodes, const auto& coarse_node) { return n_contr_nodes + coarse_node.size(); } )
								- coarse_nodes.size()
						),
						node_first_out(fine_hg.numNodes() - num_contracted_nodes + 1, 0)
		{
			if (fine_hg.hasCoarseNodes()) throw std::runtime_error("At the moment we do not support repeated contraction.");

			std::vector<nodeid> hyperedge_sizes;
			std::vector<nodeid> fine2coarse(fine_hg.numNodes(), invalid_node);
			{	//build mapping to keep hyperedge ordering and pin ordering within hyperedges.
				nodeid coarse_node = 0;
				for (const auto& sn : coarse_nodes) {
					coarse_node_weights.push_back(static_cast<nodeid>(sn.size()));
					for (nodeid x : sn)
						fine2coarse[x] = coarse_node;
					coarse_node++;
				}
				for (nodeid u = 0; u < fine_hg.numNodes(); u++)
					if (fine2coarse[u] == invalid_node)
						fine2coarse[u] = coarse_node++;
				assert(coarse_node == numNodes());
			}

			{	//write coarse pins of hyperedges into pins vector and throw away single pin hyperedges.
				boost::dynamic_bitset<> coarse_pin_inserted(numNodes());
				for (netid e = 0; e < fine_hg.numHyperedges(); e++) {
					std::size_t he_begin = pins.size();
					for (nodeid fine_pin : fine_hg.pinsOf(e)) {
						nodeid coarse_pin = fine2coarse[fine_pin];
						if (!coarse_pin_inserted[coarse_pin]) {
							coarse_pin_inserted.set(coarse_pin);
							pins.push_back(coarse_pin);
						}
						else {
							num_contracted_pins++;
						}
					}
					std::size_t num_coarse_pins_of_e = pins.size() - he_begin;
					for (nodeid coarse_pin : const_range(pins, he_begin, pins.size()))
						coarse_pin_inserted.reset(coarse_pin);
					if (num_coarse_pins_of_e == 1) {
						num_single_pin_hyperedges++;
						pins.pop_back();
					}
					if (num_coarse_pins_of_e > 1)
						hyperedge_sizes.push_back(static_cast<nodeid>(num_coarse_pins_of_e));
				}
				num_contracted_pins += num_single_pin_hyperedges;	//these were not counted above.
			}

			incident_hyperedges = std::vector<netid>(pins.size());
			hyperedge_first_out = std::vector<nodeid>(hyperedge_sizes.size() + 1, 0);
			init(hyperedge_sizes);

			assert(totalNodeWeight() == fine_hg.totalNodeWeight());
			assert(totalNodeWeight() == fine_hg.numNodes());
			assert(numNodes() == fine_hg.numNodes() - num_contracted_nodes);
			assert(numHyperedges() == fine_hg.numHyperedges() - num_single_pin_hyperedges);
			assert(totalHyperedgeWeight() == fine_hg.numHyperedges());
			assert(numPins() == fine_hg.numPins() - num_contracted_pins);
			assert(totalNumPins() == fine_hg.numPins());
			assert(numCoarseNodes() == coarse_nodes.size());
#ifndef NDEBUG
			for (nodeid u = 0; u < numNodes(); u++) {
				assert( (nodeWeight(u) == 1) == !isNodeCoarse(u) );
			}
#endif
		}

		inline nodeid numNodes() const { return static_cast<nodeid>(node_first_out.size() - 1); }
		inline netid numHyperedges() const { return static_cast<nodeid>(hyperedge_first_out.size() - 1); }
		inline std::size_t numPins() const { return pins.size(); }

		inline nodeid totalNodeWeight() const { return numNodes() + num_contracted_nodes; }
		inline netid totalHyperedgeWeight() const { return numHyperedges() + num_single_pin_hyperedges; }
		inline std::size_t totalNumPins() const { return numPins() + num_contracted_pins; }

		inline nodeid pinCount(netid e) const { return hyperedge_first_out[e+1] - hyperedge_first_out[e]; }
		inline netid degree(nodeid u) const { return node_first_out[u+1] - node_first_out[u]; }

		inline std::vector<netid>::const_iterator beginHyperedges(nodeid u) const { return incident_hyperedges.cbegin() + node_first_out[u]; }
		inline std::vector<netid>::const_iterator endHyperedges(nodeid u) const { return incident_hyperedges.cbegin() + node_first_out[u+1]; }
		inline std::vector<netid>::const_iterator endHyperedges() const { return incident_hyperedges.cend(); }

		inline std::vector<nodeid>::const_iterator beginPins(netid e) const { return pins.cbegin() + hyperedge_first_out[e]; }
		inline std::vector<nodeid>::const_iterator endPins(netid e) const { return pins.cbegin() + hyperedge_first_out[e+1]; }
		inline std::vector<nodeid>::const_iterator endPins() const { return pins.cend(); }

		const std::vector<nodeid>& getPins() const { return pins; }
		const std::vector<nodeid>& getHyperedgeFirstOut() const { return hyperedge_first_out; }

		inline const_range<netid> hyperedgesOf(nodeid u) const { return { beginHyperedges(u) , endHyperedges(u) }; }
		inline const_range<nodeid> pinsOf(netid e) const { return { beginPins(e), endPins(e) }; }

		inline bool hasCoarseNodes() const { return !coarse_node_weights.empty(); }
		inline nodeid numCoarseNodes() const { return static_cast<nodeid>(coarse_node_weights.size()); }
		inline bool isNodeCoarse(nodeid u) const { return u < numCoarseNodes(); }
		inline nodeid nodeWeight(nodeid u) const { return isNodeCoarse(u) ? coarse_node_weights[u] : 1; }

		//we need this only for hypergraphs without contractions.
		static HypergraphWithNodeWeightsWithoutHyperedgeWeights buildWithoutLargeHyperedges(nodeid num_nodes,
																				   std::vector<nodeid>& hyperedge_sizes,
																				   std::vector<nodeid> pins,
																				   nodeid max_hyperedge_size,
																				   netid& n_deleted_hyperedges)
		{
			n_deleted_hyperedges = 0;
			nodeid keep_pin_index = 0;
			nodeid global_pin_index = 0;
			netid keep_hyperedge_index = 0;
			for (netid e = 0; e < hyperedge_sizes.size(); e++) {
				if (hyperedge_sizes[e] <= max_hyperedge_size) {
					hyperedge_sizes[keep_hyperedge_index++] = hyperedge_sizes[e];
					for (nodeid i = 0; i < hyperedge_sizes[e]; i++) {
						pins[keep_pin_index++] = pins[global_pin_index++];
					}
				}
				else {
					global_pin_index += hyperedge_sizes[e];
					n_deleted_hyperedges++;
				}
			}
			hyperedge_sizes.resize(keep_hyperedge_index);
			hyperedge_sizes.shrink_to_fit();
			pins.resize(keep_pin_index);
			pins.shrink_to_fit();
			return HypergraphWithNodeWeightsWithoutHyperedgeWeights(num_nodes, hyperedge_sizes, std::move(pins));
		}
	};

	using Hypergraph = HypergraphWithNodeWeightsWithoutHyperedgeWeights;

}//namespace hyper
