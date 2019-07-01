#pragma once

#include "../definitions.h"
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>

namespace hyper {
	class Partition {
	private:
		std::vector<partitionid> part;
		netid cut = invalid_hyperedge;
		partitionid max_id_if_not_modified = invalid_block;
	public:
		explicit Partition() = default;
		explicit Partition(std::size_t size) : part(size, 0) { }
		explicit Partition(std::vector<partitionid> opart) : part(std::move(opart)) {}


		std::size_t size() { return part.size(); }

		nodeid getSetSize(partitionid p) {
			if (p > maxId()) {
				throw std::runtime_error("Invalid partition id p=" + std::to_string(p) + ".");
			}
			nodeid res = 0;
			for (partitionid px : part) { res += static_cast<nodeid>(px == p); }
			return res;
		}

		std::vector<nodeid> getSetSizes() {
			partitionid z = maxId();
			std::vector<nodeid> res(z+1, 0);
			for (partitionid b : part) { res[b]++; }
			return res;
		}

		std::vector<nodeid> membersOf(partitionid p) {
			std::vector<nodeid> res;
			for (nodeid u = 0; u < part.size(); u++) { if (part[u] == p) res.push_back(u); }
			return res;
		}

		std::pair<partitionid, std::vector<std::vector<nodeid>>> membersOfAllSubsets() {
			partitionid z = maxId();
			std::vector<std::vector<nodeid>> res(z+1);
			for (nodeid u = 0; u < part.size(); u++) { res[part[u]].push_back(u); }
			return std::make_pair(z,res);
		}

		std::vector<std::pair<partitionid, std::vector<nodeid>>> membersAndIdsOfAllSubsets() {
			std::vector<std::pair<partitionid, std::vector<nodeid>>> res(maxId() + 1);
			partitionid i = 0; for (auto& x : res) { x.first = i++; }
			for (nodeid u = 0; u < part.size(); u++) { res[part[u]].second.push_back(u); }
			return res;
		}

		partitionid& operator[](std::size_t idx) { return part[idx]; }

		const partitionid& operator[](std::size_t idx) const { return part[idx]; }

		partitionid* data() { return part.data(); }

		template<typename L> void forEntries(L l) {
			for (nodeid u = 0; u < part.size(); u++) { l(u, part[u]); }
		}

		const_range<partitionid> entries() { return { part.begin(), part.end() }; }

		void compactify() {
			partitionid fresh = 0;
			max_id_if_not_modified = 0;
			std::unordered_map<partitionid, partitionid> compactingMap;
			for (partitionid p : part) {
				auto insert_success = compactingMap.insert(std::make_pair(p, fresh)).second;
				if (insert_success) { fresh++; }
			}
			for (partitionid& p : part) {
				p = compactingMap[p];
				max_id_if_not_modified = std::max(max_id_if_not_modified, p);
			}
		}

		partitionid maxId() {
			if (part.empty()) throw std::runtime_error("No max id in empty partition");
			max_id_if_not_modified = *std::max_element(part.begin(), part.end());
			return max_id_if_not_modified;
		}

		partitionid maxIdIfNotModified() const {
			if (max_id_if_not_modified == invalid_block) throw std::runtime_error("max ID not initialized.");
			return max_id_if_not_modified;
		}

		netid computeCut(const Hypergraph& hg) {
			cut = 0;
			for (netid e = 0; e < hg.numHyperedges(); e++) {
				partitionid x = invalid_block;
				for (nodeid pin : hg.pinsOf(e)) {
					if (x == invalid_block) {
						x = part[pin];
					}
					else if (x != part[pin]) {
						cut++;
						break;
					}
				}
			}
			return cut;
		}

		netid getCut() const {
			if (cut == invalid_hyperedge) { throw std::runtime_error("Cut has not been computed yet. Call computeCut(hg)."); }
			return cut;
		}

		Partition extractRange(nodeid __begin, nodeid __end) {
			Partition sub(__end - __begin);
			std::copy(part.begin() + __begin, part.begin() + __end, sub.part.begin());
			return sub;
		}
	};

	class PartitionIntersection {
	public:
		static Partition intersect(Partition& p1, Partition& p2) {
			if (p1.size() != p2.size()) { throw std::runtime_error("p1.size() != p2.size()"); }
			partitionid z = *std::max_element(p1.entries().begin(), p1.entries().end()) + 1;
			Partition p(p1.size());
			for (nodeid u = 0; u < p1.size(); u++) { p[u] = p2[u] * z + p1[u]; }
			p.compactify();
			return p;
		}

		static Partition bulkIntersect(std::vector<Partition>& ps) {
			std::size_t remaining = ps.size();
			while (remaining > 1) {
				for (nodeid i = 0; i < remaining/2; i++) {
					ps[i] = intersect(ps[2*i], ps[2*i + 1]);
				}
				bool uneven = remaining % 2 != 0;
				remaining /= 2;
				if (uneven) {
					ps[remaining] = std::move(ps[2*remaining]);
					remaining += 1;
				}
			}
			Partition res = std::move(ps[0]);	//let's hope this works as expected
			ps.clear();
			return res;
		}
	};

	class InBlockPredicate {
	private:
		const Partition& part;
	public:
		partitionid blockid;
		InBlockPredicate(const Partition& part, const partitionid blockid) : part(part), blockid(blockid) {

		}

		bool inBlock(const nodeid x) const { return blockid == part[x]; }
	};

	class OverlapWithPartitionBlocks {
	public:
		explicit OverlapWithPartitionBlocks(Partition& p, bool use_ensemble_piercing) : use_ensemble_piercing(use_ensemble_piercing), ip(p),
																						overlaps(use_ensemble_piercing ? p.maxId() + 1 : 0,0) {}
		bool use_ensemble_piercing;
		Partition& ip;
		std::vector<nodeid> overlaps;
		void settleNode(nodeid u) { overlaps[ip[u]]++; }
		nodeid getOverlapWithNodesBlock(nodeid u) const { return overlaps[ip[u]]; }
		void swap(OverlapWithPartitionBlocks& other) {
			std::swap(overlaps, other.overlaps);
			std::swap(ip, other.ip);
		}
	};
}