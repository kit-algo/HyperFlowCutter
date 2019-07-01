#pragma once

#include "../algorithm/hla_flowcutter.h"

namespace hyper {
	class CutInfo {
	public:
		using mus_duration = std::chrono::duration<double,std::micro>;
		//CutInfo() : CutInfo(MAX_FLOW, invalid_node, invalid_node, invalid_node, -1, -1) { }
		CutInfo(flow_t c, nodeid bmin, nodeid bmax, nodeid nIso, int step = -1, int cutter_id = -1,
				mus_duration running_time_cutter = duration(0.0), mus_duration running_time_total = duration(0.0)) :
				cut(c), smallerBlock(bmin), largerBlock(bmax), nIsolatedNodes(nIso), step(step), cutter_id(cutter_id),
				running_time_cutter(running_time_cutter), running_time_total(running_time_total) {}
		flow_t cut;
		nodeid smallerBlock;
		nodeid largerBlock;
		nodeid nIsolatedNodes;
		int step;
		int cutter_id;
		std::chrono::duration<double, std::micro> running_time_cutter;
		std::chrono::duration<double, std::micro> running_time_total;

		bool isValid() { return (cutter_id != -1 && cut != MAX_FLOW); }

		static CutInfo getInfCut() {
			return CutInfo(MAX_FLOW, invalid_node, invalid_node, invalid_node);
		}

		static CutInfo getZeroCut(nodeid nNodes) {
			return CutInfo(0, 0, nNodes, 0);
		}

		static CutInfo extract(const HLAFlowCutter& flc, int cutter_id, std::chrono::duration<double, std::micro> running_time_total) {
			auto [bmin, bmax] = aux::minmax(flc.numberOfSourceReachableNodes(), flc.numberOfTargetReachableNodes());
			return CutInfo(flc.cutsize, bmin, bmax, flc.nIsolatedNodes, flc.step, cutter_id, flc.running_time, running_time_total);
		}

		static CutInfo getCutFromExternalPartitioner(flow_t cut, nodeid smallerBlock, nodeid largerBlock, std::chrono::duration<double, std::micro> running_time_total) {
			return CutInfo(cut, smallerBlock, largerBlock, 0, -2, -2, mus_duration(0.0), running_time_total);
		}
	};

	class ParetoFront {
	public:
		ParetoFront() = default;
		explicit ParetoFront(nodeid nNodes) : nNodes(nNodes), perf_balance(nNodes/2), cuts(perf_balance+1, CutInfo::getInfCut()) {
			assert(cuts.size() == perf_balance + 1);
			cuts[0] = CutInfo::getZeroCut(nNodes);
		}

		nodeid nNodes = 0;
		nodeid perf_balance = 0;
		std::vector<CutInfo> cuts;

		CutInfo getCutExactBalance(nodeid smaller_block_size) const {
			return cuts[smaller_block_size];
		}

		std::pair<nodeid,CutInfo> getCut(nodeid smaller_block_size) const {
			nodeid minind = smaller_block_size;
			for (nodeid i = smaller_block_size + 1; i <= perf_balance; i++) {
				if (cuts[i].cut <= cuts[minind].cut) { minind = i; }
			}
			return std::make_pair(minind, cuts[minind]);
		}

		std::pair<nodeid,CutInfo> getCut(double epsilon) const {
			return getCut(smallerBlockSize(epsilon));
		}

		nodeid smallerBlockSize(double epsilon) const {
			return Metrics::smallerBlockSize(nNodes, epsilon);
		}

		void writeRange(nodeid imin, nodeid imax, CutInfo& c) {
			if (imin > perf_balance) { return; }
			imax = std::min(perf_balance, imax);
			for (nodeid i = imin; i <= imax; i++) { if (c.cut < cuts[i].cut) { cuts[i] = c; } }
		}

	};

	class ParetoWrites {
	public:
		std::vector<nodeid> earliest_eligible_write;
		ParetoFront& front;	//shared

		explicit ParetoWrites(ParetoFront& front) : earliest_eligible_write(front.cuts.size()), front(front) {
			std::iota(earliest_eligible_write.begin(), earliest_eligible_write.end(), 0);
		}

		void writeRange(nodeid imin, nodeid imax, CutInfo& c) {
			if (imin > front.perf_balance) { return; }
			imax = std::min(front.perf_balance, imax);
			assert(imin < earliest_eligible_write.size());
			assert(imax < earliest_eligible_write.size());
			for (nodeid i = earliest_eligible_write[imin]; i <= imax; i = earliest_eligible_write[std::min(front.perf_balance, i+1)]) {	//perf_balance is acting as sentinel here. eew[perf_balance] is set to perf_balance + 1 after one write operation
				assert(i < front.cuts.size());
				assert(i < earliest_eligible_write.size());
				if (c.cut < front.cuts[i].cut) {
					front.cuts[i] = c;
				}
				earliest_eligible_write[i] = std::max(earliest_eligible_write[i], imax + 1);
			}
		}

		void writeAllAvailableCuts(const HLAFlowCutter& flc, int cutter_id, duration total_rt = duration(0.0)) {
			auto[bmin,bmax] = aux::minmax(flc.numberOfSourceReachableNodes(), flc.numberOfTargetReachableNodes());
			//smaller reachable set
			CutInfo c = CutInfo::extract(flc, cutter_id, total_rt);
			writeRange(bmin, bmin + flc.nIsolatedNodes, c);
			writeRange(bmax, bmax + flc.nIsolatedNodes, c);
			writeRange(bmin + flc.numberOfUnreachableNodes() - flc.nIsolatedNodes, bmin + flc.numberOfUnreachableNodes(), c);
			//bmax + unreachables - <iso> is unecessary, since if bmax + unreachables - nIso <= desired_smaller_block_size
			// --> (n - bmin) - nIso <= desired_smaller_block_size --> bmin + nIso >= n - desired_smaller_blocksize >= n/2 >= desired_smaller_blocksize
		}

		void resetEarliestEligibleWrites() {
			std::iota(earliest_eligible_write.begin(), earliest_eligible_write.end(), 0);
		}

	};
}
