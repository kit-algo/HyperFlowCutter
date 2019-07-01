//
// Created by lars on 20.02.18.
//

#pragma once

#include <cstdint>
#include <limits>
#include <vector>
#include <cassert>
#include "../definitions.h"
#include <iostream>
#include "boost/dynamic_bitset.hpp"

namespace hyper {


	template<typename id_t, typename timestamp_t=uint16_t>
	class mixed_assimilation_timestamp_set {
	public:
		explicit mixed_assimilation_timestamp_set(id_t size = 0) :
				generation(INITIAL_TIMESTAMP), timestamps(size, DEFAULT_FALSE_TIMESTAMP), _num_elements(0), _num_assimilated_elements(0),
				_num_hypernodes(0), _num_assimilated_hypernodes(0) {}

		inline size_t max_size() const {
			return timestamps.size();
		}

		inline id_t num_elements() const {
			return _num_elements;
		}

		inline id_t num_hypernodes() const {
			return _num_hypernodes;
		}

		inline id_t num_assimilated_hypernodes() const {
			return _num_assimilated_hypernodes;
		}

		inline id_t num_assimilated_elements() const {
			return _num_assimilated_elements;
		}

		inline bool contains(id_t id) const {
			assert(id < timestamps.size());
            return timestamps[id] >= generation;
		}

		inline bool is_assimilated(id_t id) const {
            assert(id < timestamps.size());
			return timestamps[id] == ASSIMILATED_TIMESTAMP;
		}

		inline void assimilate(id_t id, bool is_hypernode) {
			//reduce branching
			_num_assimilated_elements += static_cast<id_t>(!is_assimilated(id));
            _num_assimilated_hypernodes += static_cast<id_t>(!is_assimilated(id) && is_hypernode);
            _num_elements += static_cast<id_t>(!contains(id));
            _num_hypernodes += static_cast<id_t>(!contains(id) && is_hypernode);
			timestamps[id] = ASSIMILATED_TIMESTAMP;
        }

		inline void add_but_dont_assimilate(id_t id, bool is_hypernode) {
			//reduce branching
			_num_elements += static_cast<id_t>(!contains(id));
			_num_hypernodes += static_cast<id_t>(!contains(id) && is_hypernode);
			timestamps[id] = std::max(generation, timestamps[id]);
		}

		void reset_all_except_assimilated() {
 			if (generation < ASSIMILATED_TIMESTAMP - 1) {
				generation++;
			}
			else {
				generation = INITIAL_TIMESTAMP;
				for (id_t id = 0; id < timestamps.size(); id++) {
					if (timestamps[id] != ASSIMILATED_TIMESTAMP) {
						timestamps[id] = DEFAULT_FALSE_TIMESTAMP;
					}
				}
			}
			_num_elements = _num_assimilated_elements;
			_num_hypernodes = _num_assimilated_hypernodes;
		}

		timestamp_t current_timestamp() const {
			return generation;
		}

	protected:
		static constexpr timestamp_t ASSIMILATED_TIMESTAMP = std::numeric_limits<timestamp_t>::max();
		static constexpr timestamp_t DEFAULT_FALSE_TIMESTAMP = 0;
		static constexpr timestamp_t INITIAL_TIMESTAMP = 1;

		std::vector<timestamp_t> timestamps;
		timestamp_t generation;

		id_t _num_elements;
		id_t _num_assimilated_elements;

		id_t _num_hypernodes;
		id_t _num_assimilated_hypernodes;
	};

	template<typename timestamp_t, typename id_t, typename value_t, value_t initial_value, value_t assimilated_value, value_t false_value>
	class assimilation_timestamp_array : public mixed_assimilation_timestamp_set<timestamp_t, id_t> {
	public:
		explicit assimilation_timestamp_array(id_t size) :
				mixed_assimilation_timestamp_set<timestamp_t, id_t>(size), values(size, initial_value) {}

		inline const value_t& operator[](id_t id) const {
			if (is_assimilated(id)) {
				return assimilated_value;
			}
			else if (contains(id)) {
				return values[id];
			}
			else {
				return false_value;
			}
		}

		void set(id_t id, value_t& val, bool is_hypernode) {
			if (!is_assimilated(id)) {
				values[id] = val;
			}
		}

	protected:
		std::vector<value_t> values;
	};

	template<typename id_t, typename timestamp_t=uint16_t>
	class FastResetBitvector {
	public:
		static constexpr timestamp_t DEFAULT_FALSE_TIMESTAMP = 0;
		static constexpr timestamp_t INITIAL_TIMESTAMP = 1;

		explicit FastResetBitvector(id_t size = 0) : generation(INITIAL_TIMESTAMP), timestamps(size, DEFAULT_FALSE_TIMESTAMP), _num_elements(0) { }
		inline bool contains(id_t id) const { assert(id < timestamps.size()); return timestamps[id] == generation; }
		inline void remove(id_t id) { assert(id < timestamps.size()); _num_elements -= static_cast<id_t>( contains(id) ); timestamps[id] = DEFAULT_FALSE_TIMESTAMP; }
		inline void reset(id_t id) { remove(id); }
		inline void set(id_t id) { assert(id < timestamps.size()); _num_elements += static_cast<id_t>( !contains(id) ); timestamps[id] = generation; }
		inline void add(id_t id) { set(id); }
		inline id_t count() const { return _num_elements; }
		inline id_t numberOfContainedIds() const { return count(); }
		inline void reset() {
			if (generation <= std::numeric_limits<timestamp_t>::max()) {
				generation++;
			}
			else {
				generation = INITIAL_TIMESTAMP;
				for (id_t id = 0; id < timestamps.size(); id++) { timestamps[id] = DEFAULT_FALSE_TIMESTAMP; }
			}
			_num_elements = 0;
		}

	private:
		timestamp_t generation;
		std::vector<timestamp_t> timestamps;
		id_t _num_elements;
	};


	template<typename id_t, typename timestamp_t=uint16_t>
	class assimilation_timestamp_set {
	public:
		explicit assimilation_timestamp_set(id_t size = 0) :
				generation(INITIAL_TIMESTAMP), timestamps(size, DEFAULT_FALSE_TIMESTAMP), _num_elements(0), _num_assimilated_elements(0) {

		}

		void swap(assimilation_timestamp_set<id_t, timestamp_t>& other) {
			timestamps.swap(other.timestamps);
			std::swap(generation, other.generation);
			std::swap(_num_elements, other._num_elements);
			std::swap(_num_assimilated_elements, other._num_assimilated_elements);
		}

		inline size_t max_size() const {
			return timestamps.size();
		}

		inline id_t num_elements() const {
			return _num_elements;
		}

		inline id_t capacity() const { return timestamps.size(); }

		inline id_t num_assimilated_elements() const {
			return _num_assimilated_elements;
		}

		inline bool contains(id_t id) const {
			assert(id < timestamps.size());
			return timestamps[id] >= generation;
		}

		inline bool is_assimilated(id_t id) const {
			assert(id < timestamps.size());
			return timestamps[id] == ASSIMILATED_TIMESTAMP;
		}

		inline void assimilate(id_t id) {
			_num_assimilated_elements += static_cast<id_t>(!is_assimilated(id));
			_num_elements += static_cast<id_t>(!contains(id));
			timestamps[id] = ASSIMILATED_TIMESTAMP;

		}

		inline void add_but_dont_assimilate(id_t id) {
			assert(!contains(id));
			_num_elements++;
			timestamps[id] = generation;
		}

		void reset_all_except_assimilated() {
			if (generation < ASSIMILATED_TIMESTAMP - 1) {
				generation++;
			}
			else {
				generation = INITIAL_TIMESTAMP;
				for (id_t id = 0; id < timestamps.size(); id++) {
					if (!is_assimilated(id)) {
						timestamps[id] = DEFAULT_FALSE_TIMESTAMP;
					}
				}
			}
			_num_elements = _num_assimilated_elements;
		}

		timestamp_t current_timestamp() const {
			return generation;
		}
	protected:
		static constexpr timestamp_t ASSIMILATED_TIMESTAMP = std::numeric_limits<timestamp_t>::max();
		static constexpr timestamp_t DEFAULT_FALSE_TIMESTAMP = 0;
		static constexpr timestamp_t INITIAL_TIMESTAMP = 1;

		timestamp_t generation;
		std::vector<timestamp_t> timestamps;

		id_t _num_elements;
		id_t _num_assimilated_elements;
	};

	template<typename id_t, typename timestamp_t=uint16_t>
	class hyperedge_assimilation_timestamp_set {
	public:
		explicit hyperedge_assimilation_timestamp_set(id_t size = 0) :
				generation(INITIAL_TIMESTAMP), timestamps(size, DEFAULT_FALSE_TIMESTAMP) {

		}

		void swap(hyperedge_assimilation_timestamp_set<id_t, timestamp_t>& other) {
			timestamps.swap(other.timestamps);
			std::swap(generation, other.generation);
		}

		inline size_t max_size() const {
			return timestamps.size();
		}

		inline bool contains(id_t id) const {
			assert(id < timestamps.size());
			//return timestamps[id] == generation || timestamps[id] == ASSIMILATED_TIMESTAMP;
			return timestamps[id] >= generation;
		}

		inline bool is_assimilated(id_t id) const {
			assert(id < timestamps.size());
			return timestamps[id] == ASSIMILATED_TIMESTAMP;
		}

		inline void assimilate(id_t id) {
			assert(id < timestamps.size());
			timestamps[id] = ASSIMILATED_TIMESTAMP;
		}

		inline void add_but_dont_assimilate(id_t id) {
			assert(id < timestamps.size());
			assert(!contains(id));
			timestamps[id] = generation;
		}

		void reset_all_except_assimilated() {
			if (generation < ASSIMILATED_TIMESTAMP - 1) {
				generation++;
			}
			else {
				generation = INITIAL_TIMESTAMP;
				for (id_t id = 0; id < timestamps.size(); id++) {
					if (!is_assimilated(id)) {
						timestamps[id] = DEFAULT_FALSE_TIMESTAMP;
					}
				}
			}
		}

		inline id_t capacity() const { return timestamps.size(); }

		timestamp_t current_timestamp() const {
			return generation;
		}

	protected:
		static constexpr timestamp_t ASSIMILATED_TIMESTAMP = std::numeric_limits<timestamp_t>::max();
		static constexpr timestamp_t DEFAULT_FALSE_TIMESTAMP = 0;
		static constexpr timestamp_t INITIAL_TIMESTAMP = 1;

		timestamp_t generation;
		std::vector<timestamp_t> timestamps;
	};

	template<typename id_t>
	class BitVectorAssimilationSet {
	public:
		explicit BitVectorAssimilationSet(id_t size) : _is_assimilated(size), is_reachable(size) { }
		void swap(BitVectorAssimilationSet& other) { _is_assimilated.swap(other._is_assimilated); is_reachable.swap(other.is_reachable); }
		void reset_all_except_assimilated() { is_reachable = _is_assimilated; }
		inline void assimilate(id_t id) { _is_assimilated.set(id); is_reachable.set(id); }
		inline bool is_assimilated(id_t id) const { return _is_assimilated[id]; }
		inline void add_but_dont_assimilate(id_t id) { is_reachable.set(id); }
		inline bool contains(id_t id) const { return is_reachable[id]; }
		inline id_t num_assimilated_elements() const { return static_cast<id_t>(_is_assimilated.count()); }
		inline id_t num_elements() const { return static_cast<id_t>(is_reachable.count()); }
		inline id_t capacity() const { return static_cast<id_t>(is_reachable.size()); }
	protected:
		boost::dynamic_bitset<> _is_assimilated;
		boost::dynamic_bitset<> is_reachable;
	};

}//namespace hyper
