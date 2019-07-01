#pragma once

#include "boost/dynamic_bitset.hpp"


template<typename id_t>
class SparseResetBitvector {
private:
	std::vector<id_t> contained_ids;
	boost::dynamic_bitset<> is_contained;

	inline void addImpl(id_t t) {
		contained_ids.push_back(t);
		is_contained.set(t);
	}

public:

	explicit SparseResetBitvector(std::size_t num_elements) : is_contained(num_elements) { }

	inline bool contains(id_t t) { return is_contained[t]; }

	void reset() {
		if (contained_ids.size() > is_contained.size()/boost::dynamic_bitset<>::bits_per_block) {
			is_contained.reset();
		}
		else {
			for (id_t t : contained_ids)
				is_contained.reset(t);
		}
		contained_ids.clear();
	}

	template<bool check_if_contained=false>
	inline void add(id_t t) {
		if (check_if_contained) {
			if (!contains(t)) {
				addImpl(t);
			}
		}
		else {
			addImpl(t);
		}
	}

	std::vector<id_t>& get_contained_ids() { return contained_ids; }
	std::size_t numberOfContainedIds() { return contained_ids.size(); }
	std::size_t size() { return contained_ids.size(); }

	std::vector<id_t> extractIds() { return std::move(contained_ids); }
	boost::dynamic_bitset<> extractBitset() { return std::move(is_contained); }
};

template<typename id_t, typename count_t>
class SparseResetCountVector {
private:
	std::vector<id_t> contained_ids;
	std::vector<count_t> values;
public:
	explicit SparseResetCountVector(std::size_t num_elements) : values(num_elements, 0) { }
	inline bool contains(id_t t) { return values[t] > 0; }
	inline count_t& operator[](id_t t) { return values[t]; }
	void reset() {
		for (id_t t : contained_ids)
			values[t] = 0;
		contained_ids.clear();
	}
	inline void add(id_t t, count_t val) {
		if (values[t] == 0)
			contained_ids.push_back(t);
		values[t] += val;
	}
	std::vector<id_t>& get_contained_ids() { return contained_ids; }
	id_t numberOfContainedIds() { return static_cast<id_t>(contained_ids.size()); }

	void sortByKey() { std::sort(contained_ids.begin(), contained_ids.end()); }

	template<typename func_handle_kv_pair>
	void forKeyValuePairs(const func_handle_kv_pair& handle_kv_pair) {
		for (id_t id : contained_ids) {
			handle_kv_pair(id, values[id]);
		}
	}
};