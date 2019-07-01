#pragma once

#include "../../datastructure/partition.h"
#include "../../algorithm/hla_flowcutter.h"
#include "../../io/hmetisreader.h"
#include "../../io/patohreader.h"

#include "../../extern/patoh_interface.h"


namespace hyper {
	class EnsemblePartitioning {
	public:
		EnsemblePartitioning() = default;
		EnsemblePartitioning(std::vector<Partition>& ps, const Hypergraph& hg) : nPartitions(ps.size()), intersectionPartition(0),
																					  cutFrequency(hg.numHyperedges(), 0) {
			if (nPartitions == 0) { return; }
			if (nPartitions > MAX_NUM_PATOH_PARTITIONS) {
				std::string error_message = "At most " + std::to_string(MAX_NUM_PATOH_PARTITIONS) + " permitted."
											+ "If you need more, adapt the PATOH_SCORE_RANGE parameter at own risk. (Make sure the scores don't collide with other scorers.)";
				throw std::runtime_error(error_message);
			}
			for (netid e = 0; e < hg.numHyperedges(); e++) {
				for (Partition& p : ps) {
					bool s0 = false,s1 = false;
					for (nodeid pin : hg.pinsOf(e)) {
						s0 |= p[pin] == 0;
						s1 |= p[pin] == 1;
						if (s0 && s1) { break; }
					}
					if (s0 && s1) {
						cutFrequency[e]++;
					}
				}
			}
			intersectionPartition = PartitionIntersection::bulkIntersect(ps);
			maxPartitionID = *std::max_element(intersectionPartition.entries().begin(), intersectionPartition.entries().end());
		}

		std::size_t nPartitions = 0;	//recommend 10
		Partition intersectionPartition;
		std::vector<uint32_t> cutFrequency;
		partitionid maxPartitionID = 0;
		std::size_t numUsedEnsembleSTPairs = 0;

		std::vector<nodeid> getOverlapWithIntersectionPartitionBlocks(std::vector<nodeid>& nodes) {
			std::vector<nodeid> overlap(maxPartitionID+1, 0);
			for (nodeid u : nodes) { overlap[intersectionPartition[u]] ++; }
			return overlap;
		}

		static EnsemblePartitioning buildPatohEnsemblePartitioning(const Hypergraph& hg, State& state) {
			if (state.piercerOptions.numEnsemblePartitions == 0) {
				return EnsemblePartitioning();
			}
			std::vector<int> seeds = Random::randomSequence<int>(state.piercerOptions.numEnsemblePartitions);
			std::vector<Partition> partitions = PaToHInterface::bisectWithPatoh(hg, seeds);
			return EnsemblePartitioning(partitions, hg);
		}

		static EnsemblePartitioning readPatohEnsemblePartitioning(const Hypergraph& hg, std::vector<std::string>& partfiles) {
			std::vector<Partition> parts;
			for (std::string& partfile : partfiles) { parts.push_back(PatohPartitionReader::read(partfile)); }
			return EnsemblePartitioning(parts, hg);
		}

		static EnsemblePartitioning readPatohEnsemblePartitioning(std::vector<std::string>& partfiles, std::string& hgfile) {
			Hypergraph hg(HMetisReader::read(hgfile));
			return readPatohEnsemblePartitioning(hg, partfiles);
		}

		static constexpr uint64_t PATOH_HYPEREDGE_MAXSCORE = uint64_t(1) << 63U;
		static constexpr uint32_t PATOH_SCORE_RANGE = 6; //~2^6 possible scores for PaToH
		static constexpr uint64_t MAX_NUM_PATOH_PARTITIONS = (uint64_t(1) << PATOH_SCORE_RANGE) - 1;
		static constexpr uint64_t PATOH_HYPEREDGE_MINSCORE = PATOH_HYPEREDGE_MAXSCORE >> PATOH_SCORE_RANGE;
		static constexpr uint64_t PATOH_NODE_MAXSCORE = uint64_t(1) << 63U;
		static constexpr uint64_t PATOH_NODE_MINSCORE = uint64_t(1) << 56U;

		inline uint64_t scoreHyperedge(netid e, const HLAFlowCutter& flc, const Hypergraph& hg) {
			return std::min(PATOH_HYPEREDGE_MAXSCORE - cutFrequency[e], PATOH_HYPEREDGE_MINSCORE);
		}

		inline uint64_t scoreNode(nodeid u, const HLAFlowCutter& flc, const Hypergraph& hg) {
			//overlap with flc's current source side.
			nodeid overlap = flc.source_side_overlap_with_ensemble_partition.getOverlapWithNodesBlock(u);
			return std::max(PATOH_NODE_MINSCORE, PATOH_NODE_MAXSCORE - overlap);
		}
	};
}