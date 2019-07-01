#pragma once
#include <vector>

namespace hyper {
	template<typename id_t>
	class PiercingQueue {
	private:
		std::vector<id_t> queue;
		id_t numberOfPiercingNodes;

		id_t qFront;
		id_t qEnd;
	public:
		explicit PiercingQueue(const id_t num_entries) : queue(num_entries) {}
		inline bool empty() { return qFront >= qEnd; }
		inline void push(id_t x) { assert(qEnd < queue.size()); queue[qEnd++] = x; }
		inline id_t pop() { assert(qFront < qEnd); return queue[qFront++]; }
		inline void reinitializeQueue() { qFront = 0; qEnd = numberOfPiercingNodes; }
		inline void insertPiercingNode(id_t x) { assert(numberOfPiercingNodes < queue.size()); queue[numberOfPiercingNodes++] = x; }
		inline void clearPiercingNodes() { numberOfPiercingNodes = 0; }
	};

}
