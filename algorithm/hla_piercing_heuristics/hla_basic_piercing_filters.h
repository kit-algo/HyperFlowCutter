#pragma once

#include <cstdint>
#include <random>
#include "../../datastructure/sparse_reset_vector.h"
#include "../hla_flowcutter.h"
#include "hla_basic_scorers.h"
#include "patoh_ensemble_partitioning.h"

namespace hyper {

/*
	class AllScorer {
	public:
		State& state;
		PiercingScorers::ThreeHopTarget tht;
		PiercingScorers::RandomTiebreak rt;
		EnsemblePartitioning& ep;

		explicit AllScorer(netid numHyperedges, State& state, EnsemblePartitioning& ep) : ep(ep),
																						  tht(state.piercerOptions.use_threehop ? numHyperedges : 0),
																						  state(state) {}

		inline uint64_t scoreHyperedge(netid piercing_hyperedge, const HLAFlowCutter& flc, const Hypergraph& hg) {
			uint64_t score = 0;
			if (state.piercerOptions.use_threehop) { score += tht.score(piercing_hyperedge, hg, flc); }
			if (state.piercerOptions.use_random) { score += rt.score(); }
			if (state.piercerOptions.use_ensemble) { score += ep.scoreHyperedge(piercing_hyperedge, flc, hg); }
			return score;
		}

		static constexpr uint64_t TWOHOP_MAXSCORE = uint64_t(1) << 55U;
		static constexpr uint64_t TWOHOP_MINSCORE = uint64_t(1) << 25U;
		inline uint64_t scoreNode(nodeid p, const HLAFlowCutter& flc, const Hypergraph& hg) {
			uint64_t score = 0;
			if (state.piercerOptions.use_threehop) {
				netid num_uncut_hyperedges_with_target_pins = flc.number_of_hyperedges_with_target_pins[p] - flc.number_of_mixed_incident_hyperedges[p];
				score += std::max(TWOHOP_MAXSCORE - num_uncut_hyperedges_with_target_pins, TWOHOP_MINSCORE);
			}
			if (state.piercerOptions.use_random) { score += rt.score(); }
			if (state.piercerOptions.use_ensemble) { score += ep.scoreNode(p, flc, hg); }
			return score;
		}
	};

	namespace PiercingFilters {
		class PiercingPipeline {
		private:
			State& state;
			AllScorer scorer;

			inline bool hyperedgeAvoidsAugmentingPath(const HLAFlowCutter &flc, const Hypergraph &hg, netid e) {
				return !flc.isCutHyperedgeReachableFromTarget(e, hg);
			}
			inline bool nodeAvoidsAugmentingPath(const HLAFlowCutter& flc, const Hypergraph& hg, nodeid u) {
				return !flc.isReachableFromTarget(u);
			}
			inline nodeid numberOfUnassimilatedPins(const HLAFlowCutter& flc, const Hypergraph& hg, netid e) {
				return flc.activePinCount(e);
			}
			inline bool doesNotBreakBlocksize(const HLAFlowCutter& flc, nodeid numberOfUnassimilatedPins) {
				return flc.numberOfSourceAssimilatedNodes() + numberOfUnassimilatedPins <= flc.totalNodeWeight()/2;
			}

		public:
			explicit PiercingPipeline(const netid numHyperedges, State& state, EnsemblePartitioning& ep) : scorer(totalHyperedgeWeight, state, ep), state(state) {}

			std::pair<bool,bool> findAvoidAugmentingPathPiercingHyperedge(const HLAFlowCutter& flc, const Hypergraph& hg, netid& piercing_hyperedge) {
				if (!state.piercerOptions.avoid_augmenting_paths) { bool res = findAnyPiercingHyperedge(flc, hg, piercing_hyperedge); return std::make_pair(res,res); }
				piercing_hyperedge = invalid_hyperedge;
				bool eligible_candidate = false;
				uint64_t max_score = 0;
				for (netid e : flc.source_cut_front) {
					nodeid unassimPins = numberOfUnassimilatedPins(flc, hg, e);
					if (doesNotBreakBlocksize(flc, unassimPins)) {
						eligible_candidate = true;
						if (hyperedgeAvoidsAugmentingPath(flc, hg, e)) {
							uint64_t score_e = scorer.scoreHyperedge(e, flc, hg);
							if (score_e > max_score) { max_score = score_e; piercing_hyperedge = e; }
						}
					}
				}
				bool found = piercing_hyperedge != invalid_hyperedge;
				return std::make_pair(eligible_candidate, found);
			}

			bool avoidAugmentingPaths() { return state.piercerOptions.avoid_augmenting_paths; }

			bool findAnyPiercingHyperedge(const HLAFlowCutter& flc, const Hypergraph& hg, netid& piercing_hyperedge) {
				piercing_hyperedge = invalid_hyperedge;
				uint64_t max_score = 0;
				for (netid e : flc.source_cut_front) {
					nodeid unassimPins = numberOfUnassimilatedPins(flc, hg, e);
					if (doesNotBreakBlocksize(flc, unassimPins)) {
						uint64_t score_e = scorer.scoreHyperedge(e, flc, hg);
						if (score_e > max_score) { max_score = score_e; piercing_hyperedge = e; }
					}
				}
				return piercing_hyperedge != invalid_hyperedge;
			}

			nodeid findAvoidAugmentingPathPiercingNode(const HLAFlowCutter& flc, const Hypergraph& hg) {
				if (!state.piercerOptions.avoid_augmenting_paths) { return findAnyPiercingNode(flc, hg); }
				uint64_t max_score = 0;
				nodeid piercing_node = invalid_node;
				for (nodeid p : flc.source_border_vertices) {
					if (nodeAvoidsAugmentingPath(flc, hg, p)) {
						uint64_t score = scorer.scoreNode(p, flc, hg);
						if (score > max_score) { max_score = score; piercing_node = p; }
					}
				}
				return piercing_node;
			}

			nodeid findAnyPiercingNode(const HLAFlowCutter& flc, const Hypergraph& hg) {
				uint64_t max_score = 0;
				nodeid piercing_node = invalid_node;
				for (nodeid p : flc.source_border_vertices) {
					uint64_t score = scorer.scoreNode(p, flc, hg);
					if (score > max_score) { max_score = score; piercing_node = p; }
				}
				return piercing_node;
			}
		};

	};
*/
	//Too lazy to make properly generic.
	namespace IgnoreOptionsPiercing {

	class IgnoreOptionsScorer {
	public:
		EnsemblePartitioning& ep;
		PiercingScorers::ThreeHopTarget tht;
		PiercingScorers::RandomTiebreak rt;

		explicit IgnoreOptionsScorer(netid numHyperedges, EnsemblePartitioning& ep) : ep(ep), tht(numHyperedges) {}

		inline uint64_t scoreHyperedge(netid piercing_hyperedge, const HLAFlowCutter& flc, const Hypergraph& hg) {
			return rt.score();
			//return tht.score(piercing_hyperedge, hg, flc) +rt.score() + ep.scoreHyperedge(piercing_hyperedge, flc, hg);
		}

		static constexpr uint64_t TWOHOP_MAXSCORE = uint64_t(1) << 55U;
		static constexpr uint64_t TWOHOP_MINSCORE = uint64_t(1) << 25U;
		inline uint64_t scoreNode(nodeid p, const HLAFlowCutter& flc, const Hypergraph& hg) {
			return rt.score();
			//netid num_uncut_hyperedges_with_target_pins = flc.number_of_hyperedges_with_target_pins[p] - flc.number_of_mixed_incident_hyperedges[p];
			//return std::max(TWOHOP_MAXSCORE - num_uncut_hyperedges_with_target_pins, TWOHOP_MINSCORE) + rt.score();// + ep.scoreNode(p, flc, hg);
		}
	};

		class PiercingPipeline {
		private:
			IgnoreOptionsScorer scorer;

			inline bool hyperedgeAvoidsAugmentingPath(const HLAFlowCutter &flc, netid e) {
				return !flc.isCutHyperedgeReachableFromTarget(e);
			}
			inline bool nodeAvoidsAugmentingPath(const HLAFlowCutter& flc, nodeid u) {
				return !flc.isReachableFromTarget(u);
			}
			inline nodeid numberOfUnassimilatedPins(const HLAFlowCutter& flc, netid e) {
				return flc.activePinCount(e);
			}
			inline bool doesNotBreakBlocksize(const HLAFlowCutter& flc, nodeid numberOfUnassimilatedPins) {
				return flc.numberOfSourceAssimilatedNodes() + numberOfUnassimilatedPins <= flc.totalNodeWeight()/2;
			}

		public:
			explicit PiercingPipeline(const netid numHyperedges, State& state, EnsemblePartitioning& ep) : scorer(numHyperedges, ep) {}
			std::pair<bool,bool> findAvoidAugmentingPathPiercingHyperedge(const HLAFlowCutter& flc, const Hypergraph& hg, netid& piercing_hyperedge) {
				piercing_hyperedge = invalid_hyperedge;
				bool eligible_candidate = false;
				uint64_t max_score = 0;
				for (netid e : flc.source_cut_front) {
					nodeid unassimPins = numberOfUnassimilatedPins(flc, e);
					if (doesNotBreakBlocksize(flc, unassimPins)) {
						eligible_candidate = true;
						if (hyperedgeAvoidsAugmentingPath(flc, e)) {
							uint64_t score_e = scorer.scoreHyperedge(e, flc, hg);
							if (score_e > max_score) { max_score = score_e; piercing_hyperedge = e; }
						}
					}
				}
				bool found = piercing_hyperedge != invalid_hyperedge;
				return std::make_pair(eligible_candidate, found);
			}

			bool findAnyPiercingHyperedge(const HLAFlowCutter& flc, const Hypergraph& hg, netid& piercing_hyperedge) {
				piercing_hyperedge = invalid_hyperedge;
				uint64_t max_score = 0;
				for (netid e : flc.source_cut_front) {
					nodeid unassimPins = numberOfUnassimilatedPins(flc, e);
					if (doesNotBreakBlocksize(flc, unassimPins)) {
						uint64_t score_e = scorer.scoreHyperedge(e, flc, hg);
						if (score_e > max_score) { max_score = score_e; piercing_hyperedge = e; }
					}
				}
				return piercing_hyperedge != invalid_hyperedge;
			}

			nodeid findAvoidAugmentingPathPiercingNode(const HLAFlowCutter& flc, const Hypergraph& hg) {
				uint64_t max_score = 0;
				nodeid piercing_node = invalid_node;
				for (nodeid p : flc.source_border_vertices) {
					if (nodeAvoidsAugmentingPath(flc, p)) {
						uint64_t score = scorer.scoreNode(p, flc, hg);
						if (score > max_score) { max_score = score; piercing_node = p; }
					}
				}
				return piercing_node;
			}

			nodeid findAnyPiercingNode(const HLAFlowCutter& flc, const Hypergraph& hg) {
				nodeid piercing_node = invalid_node;
				if (!flc.source_border_vertices.empty()) {
					piercing_node = Random::sampleOne(flc.source_border_vertices);
				}
				/*
				uint64_t max_score = 0;
				for (nodeid p : flc.source_border_vertices) {
					uint64_t score = scorer.scoreNode(p, flc, hg);
					if (score > max_score) { max_score = score; piercing_node = p; }
				}
				 */
				return piercing_node;
			}
	};

	}
}