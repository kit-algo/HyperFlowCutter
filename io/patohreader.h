#pragma once

#include <vector>
#include "../definitions.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "../datastructure/partition.h"

namespace hyper {
	class PatohPartitionReader {
	public:
		static Partition read(std::string& file) {
			std::vector<partitionid> part;
			std::ifstream f(file);
			if (!f) { throw std::runtime_error("File " + file + " not found."); }
			std::string line;
			while (std::getline(f, line)) {
				std::istringstream iss(line);
				partitionid p;
				while (iss >> p) {
					part.push_back(p);
				}
			}
			f.close();
			return Partition(part);
		}
	};


}