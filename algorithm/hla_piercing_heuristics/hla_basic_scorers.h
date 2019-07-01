#pragma once

#include <cstdint>
#include <random>
#include "../../datastructure/sparse_reset_vector.h"
#include "../hla_flowcutter.h"

namespace hyper {
namespace PiercingScorers {

	class RandomTiebreak {
	public:
		static constexpr uint64_t RECOMMENDED_MIN_SCORE = uint64_t(1);
		static constexpr uint64_t RECOMMENDED_MAX_SCORE = (uint64_t(1) << 23U) - 1;
		uint64_t MIN_SCORE = RECOMMENDED_MIN_SCORE;
		uint64_t MAX_SCORE = RECOMMENDED_MAX_SCORE;
		explicit RandomTiebreak() :  distribution(MIN_SCORE, MAX_SCORE) {}
		void modify_range(uint64_t min_score, uint64_t max_score) { MIN_SCORE = min_score; MAX_SCORE = max_score; distribution = std::uniform_int_distribution<uint64_t>(MIN_SCORE, MAX_SCORE); }
		inline uint64_t score() { return distribution(Random::getRNG()); }
		inline uint64_t score(netid e, const Hypergraph& hg, const HLAFlowCutter& flc) { return score(); }
	private:
		std::uniform_int_distribution<uint64_t> distribution;
	};

	class ThreeHopTarget {
	private:

		//typedef SparseResetBitvector<netid> FastResetHyperedgeSet;
		typedef FastResetBitvector<netid, uint8_t> FastResetHyperedgeSet;
		FastResetHyperedgeSet unique_adjacent_hyperedges;
		std::vector<nodeid> high_degree_pins;
	public:
		netid maxSampleDegree = 2500;
		//netid maxSampleDegree = 250;
		nodeid maxSamplePincount = 4000;
		//nodeid maxSamplePincount = 400;
		explicit ThreeHopTarget(const netid numHyperedges) : unique_adjacent_hyperedges(numHyperedges) {}
		inline netid numberOfIncidentHyperedgesWithTargetPins(netid piercing_hyperedge, const Hypergraph& hg, const HLAFlowCutter& flc) {
			high_degree_pins.clear();
			unique_adjacent_hyperedges.reset();
			nodeid psample = 0;
			for (nodeid pin : hg.pinsOf(piercing_hyperedge)) {
				if (flc.isSource(pin)) continue;
				if (psample++ > maxSamplePincount) { break; }
				if (hg.degree(pin) > maxSampleDegree) {
					high_degree_pins.push_back(pin);
				}
				else {
					for (netid e : hg.hyperedgesOf(pin)) {
						if (!unique_adjacent_hyperedges.contains(e) && e != piercing_hyperedge && flc.hasTargetAssimilatedPin(e) && !flc.hasSourceAssimilatedPin(e)) {
							unique_adjacent_hyperedges.add(e);
						}
					}
				}
			}
			auto result = static_cast<netid>(unique_adjacent_hyperedges.numberOfContainedIds());
			for (nodeid pin : high_degree_pins) {
				netid sample = 0;
				netid additional_per_high_degree_pin = 0;
				for (netid e : hg.hyperedgesOf(pin)) {
					if (!unique_adjacent_hyperedges.contains(e) && e != piercing_hyperedge && flc.hasTargetAssimilatedPin(e) && !flc.hasSourceAssimilatedPin(e)) {
						unique_adjacent_hyperedges.add(e);
						additional_per_high_degree_pin++;
					}
					if (++sample > maxSampleDegree) { break; }
				}
				result += additional_per_high_degree_pin * hg.degree(pin) / maxSampleDegree;
			}
			if (flc.activePinCount(piercing_hyperedge) > maxSamplePincount) return result * flc.activePinCount(piercing_hyperedge) / maxSamplePincount;
			else return result;
		}
		static constexpr uint64_t MAX_SCORE = uint64_t(1) << 55U;
		static constexpr uint64_t MIN_SCORE = uint64_t(1) << 25U;
		inline uint64_t score(netid e, const Hypergraph& hg, const HLAFlowCutter& flc) { return std::max(MAX_SCORE - numberOfIncidentHyperedgesWithTargetPins(e, hg, flc), MIN_SCORE); }
	};

}
}