#pragma once

#include <vector>
#include "../definitions.h"

template<typename id_t>
class LayeredQueue {
private:
	std::vector<id_t> queue;
public:
	id_t layerfront, layerend, qfront, qend;
	explicit LayeredQueue(const id_t num_elements) : queue(num_elements), layerfront(0), layerend(0), qfront(0), qend(0) {}
	inline void clear() { layerfront = 0; layerend = 0; qfront = 0; qend = 0; }
	inline bool empty() const { return qfront == qend; }
	inline bool currentLayerEmpty() const { return qfront == layerend; }
	inline id_t pop() { return queue[qfront++]; }
	inline id_t previousLayerPop() { return queue[layerfront++]; }
	inline void finishNextLayer() { layerend = qend; }
	inline void push(const id_t x) { assert(qend < queue.size()); queue[qend++] = x; }
	inline bool previousLayerEmpty() const { return layerfront == layerend; }
	inline id_t capacity() const { return static_cast<id_t>(queue.size()); }
	inline std::vector<id_t>& data() { return queue; }
	template<typename Func> inline void forAllEverContainedElements(Func f) { for (id_t i = 0; i < qend; i++) { f(queue[i]); } }
	inline const_range<id_t> range(id_t __begin, id_t __end) { return { queue.begin() + __begin, queue.cbegin() + __end }; }
	inline const_range<id_t> currentLayer() { return range(layerfront, qend); }
	inline id_t numberOfPushedElements() const { return qend; }

	inline id_t elementAt(const id_t t) const { return queue[t]; }
	inline void setTo(const id_t pos, id_t element) { queue[pos] = element; }

	inline id_t swapFrontToPositionAndPop(id_t pos) {
		std::swap(queue[pos], queue[qfront]);
		return pop();
	}

	std::vector<id_t> extract() { return std::move(queue); }

	template<bool resize = true>
	void inject(std::vector<id_t> external_queue, id_t num_elements) {
		queue = std::move(external_queue);
		layerfront = 0;
		layerend = 0;
		qfront = 0;
		qend = queue.size();
		if (resize)
			queue.resize(num_elements);
	}
};

namespace hyper {
	class LayeredBFS {
	public:
		static void run(const Hypergraph& hg, nodeid s, LayeredQueue<nodeid>& queue, boost::dynamic_bitset<>& node_visited, boost::dynamic_bitset<>& he_used, std::vector<nodeid>& layer_ranges) {
			queue.clear(); node_visited.reset(); he_used.reset(); layer_ranges.clear();
			queue.push(s); layer_ranges.push_back(0);
			node_visited.set(s);
			while (!queue.empty()) {
				queue.finishNextLayer();
				layer_ranges.push_back(queue.layerend);
				while (!queue.currentLayerEmpty()) {
					nodeid u = queue.pop();
					for (netid e : hg.hyperedgesOf(u)) {
						if (!he_used[e]) {
							he_used.set(e);
							for (nodeid v : hg.pinsOf(e)) {
								if (!node_visited[v]) {
									queue.push(v); node_visited.set(v);
								}
							}
						}
					}
				}
			}
		}
	};
}