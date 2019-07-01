#pragma once

#include <vector>
#include <cassert>

template<typename index_t, typename element_t>
class FixedCapacityStack {
private:
	std::vector<element_t> stack;
	index_t __size;
public:
	explicit FixedCapacityStack(const index_t num_elements) : stack(num_elements), __size(0) {}
	inline void clear() { __size = 0; }
	inline bool empty() const { return __size == 0; }
	inline element_t pop() { assert(!empty()); return stack[--__size]; }
	inline element_t top() { assert(!empty()); return stack[__size - 1]; }
	inline void push(const element_t& x) { assert(__size < stack.size()); stack[__size++] = x; }
	inline element_t elementAt(const index_t t) const { return stack[t]; }
	inline index_t size() { return __size; }

	inline index_t capacity() const { return static_cast<index_t>(stack.size()); }
	inline std::vector<element_t>& data() { return stack; }
};
