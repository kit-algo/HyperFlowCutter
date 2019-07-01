#ifndef ID_QUEUE_H
#define ID_QUEUE_H

#include "../definitions.h"

#include <algorithm>
#include <cassert>
#include <vector>

namespace hyper {
	struct IDKeyPair {
		nodeid id;
		flow_t key;
	};

//! A priority queue where the elements are IDs from 0 to id_count-1 where id_count is a number that
//! is set in the constructor.
//! The elements are sorted by integer keys.
	class MinIDQueue {
	private:
		static const nodeid tree_arity = 4;

	public:
		MinIDQueue() : heap_size(0) {}

		explicit MinIDQueue(nodeid id_count)
				: id_pos(id_count, invalid_node), heap(id_count), heap_size(0) {}

		//! Returns whether the queue is empty. Equivalent to checking whether size() returns 0.
		bool empty() const { return heap_size == 0; }

		//! Returns the number of elements in the queue.
		nodeid size() const { return heap_size; }

		//! Returns the id_count value passed to the constructor.
		nodeid id_count() const { return static_cast<nodeid>(id_pos.size()); }

		//! Checks whether an element is in the queue.
		bool contains_id(nodeid id) {
			assert(id < id_count());
			return id_pos[id] != invalid_node;
		}

		//! Removes all elements from the queue.
		void clear() {
			for (nodeid i = 0; i < heap_size; ++i)
				id_pos[heap[i].id] = invalid_node;
			heap_size = 0;
		}

		friend void swap(MinIDQueue &l, MinIDQueue &r) {
			using std::swap;
			swap(l.id_pos, r.id_pos);
			swap(l.heap, r.heap);
			swap(l.heap, r.heap);
			swap(l.heap_size, r.heap_size);
		}

		//! Returns the current key of an element.
		//! Undefined if the element is not part of the queue.
		flow_t get_key(nodeid id) const {
			assert(id < id_count());
			assert(id_pos[id] != invalid_node);
			return heap[id_pos[id]].key;
		}

		//! Returns the smallest element key pair without removing it from the queue.
		IDKeyPair peek() const {
			assert(!empty());
			return heap.front();
		}

		//! Returns the smallest element key pair and removes it form the queue.
		IDKeyPair pop() {
			assert(!empty());
			--heap_size;
			std::swap(heap[0].key, heap[heap_size].key);
			std::swap(heap[0].id, heap[heap_size].id);
			id_pos[heap[0].id] = 0;
			id_pos[heap[heap_size].id] = invalid_node;

			move_down_in_tree(0);
			return heap[heap_size];
		}

		//! Inserts a element key pair.
		//! Undefined if the element is part of the queue.
		void push(IDKeyPair p) {
			assert(p.id < id_count());
			assert(!contains_id(p.id));

			nodeid pos = heap_size;
			++heap_size;
			heap[pos] = p;
			id_pos[p.id] = pos;
			move_up_in_tree(pos);
		}

		//! Updates the key of an element if the new key is smaller than the old key.
		//! Does nothing if the new key is larger.
		//! Undefined if the element is not part of the queue.
		bool decrease_key(IDKeyPair p) {
			assert(p.id < id_count());
			assert(contains_id(p.id));

			nodeid pos = id_pos[p.id];

			if (heap[pos].key > p.key) {
				heap[pos].key = p.key;
				move_up_in_tree(pos);
				return true;
			} else {
				return false;
			}
		}

		//! Updates the key of an element if the new key is larger than the old key.
		//! Does nothing if the new key is smaller.
		//! Undefined if the element is not part of the queue.
		bool increase_key(IDKeyPair p) {
			assert(p.id < id_count());
			assert(contains_id(p.id));

			nodeid pos = id_pos[p.id];

			if (heap[pos].key < p.key) {
				heap[pos].key = p.key;
				move_down_in_tree(pos);
				return true;
			} else {
				return false;
			}
		}

	private:
		void move_up_in_tree(nodeid pos) {
			while (pos != 0) {
				nodeid parent = (pos - 1) / tree_arity;
				if (heap[parent].key > heap[pos].key) {
					std::swap(heap[pos], heap[parent]);
					std::swap(id_pos[heap[pos].id], id_pos[heap[parent].id]);
				}
				pos = parent;
			}
		}

		void move_down_in_tree(nodeid pos) {
			for (;;) {
				nodeid first_child = tree_arity * pos + 1;
				if (first_child >= heap_size)
					return; // no children
				nodeid smallest_child = first_child;
				for (nodeid c = first_child + 1;
					 c < std::min(tree_arity * pos + tree_arity + 1, heap_size); ++c) {
					if (heap[smallest_child].key > heap[c].key) {
						smallest_child = c;
					}
				}

				if (heap[smallest_child].key >= heap[pos].key)
					return; // no child is smaller

				std::swap(heap[pos], heap[smallest_child]);
				std::swap(id_pos[heap[pos].id], id_pos[heap[smallest_child].id]);
				pos = smallest_child;
			}
		}

		std::vector<nodeid> id_pos;
		std::vector<IDKeyPair> heap;

		nodeid heap_size;
	};
}//namespace
#endif
