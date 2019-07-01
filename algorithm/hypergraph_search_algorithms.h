#pragma once

namespace hyper {
	class HLAFlowCutter;

	class SearchAlgorithm {
	protected:
		const Hypergraph& hg;
		explicit SearchAlgorithm(const Hypergraph& hg, bool build_datastructures_during_grow_reachable) : hg(hg) { }
	public:
		class CutterSpecificDatastructures {
		public:
			virtual ~CutterSpecificDatastructures() = default;
			auto clone() const { return std::unique_ptr<CutterSpecificDatastructures>(clone_impl()); }
		protected:
				virtual CutterSpecificDatastructures* clone_impl() const = 0;
		};

		virtual ~SearchAlgorithm() = default;
		virtual std::unique_ptr<CutterSpecificDatastructures> obtainOwnedCutterSpecificDatastructures() const = 0;
		virtual void growAssimilated(HLAFlowCutter& flc) = 0;
		virtual void growReachable(HLAFlowCutter& flc) = 0;
		virtual void growReachableWithoutReset(HLAFlowCutter& flc) = 0;
		virtual flow_t exhaustFlow(HLAFlowCutter& flc) = 0;
	};

}//namespace
