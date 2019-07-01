#include "hmetisreader.h"

void hyper::read_hmetis(std::string& filename, hypergraph_info& hg) {
	std::ifstream f(filename);
	if (!f) { throw std::runtime_error("File: " + filename + " not found."); }

	std::string line;
	mgetline(f,line);
	{
		std::istringstream iss(line);
		uint8_t type;
		iss >> hg.num_nets >> hg.num_nodes;
		if (!(iss >> type)) {
			type = static_cast<uint8_t>(HYPERGRAPH_TYPE::UNWEIGHTED);
		}
		hg.type = static_cast<HYPERGRAPH_TYPE >(type);
	}
	hg.pins.clear();
	hg.netsizes.clear();
	hg.num_oversized_nets = 0;

	bool edge_weights = hg.type == HYPERGRAPH_TYPE::EDGEANDNODEWEIGHTS || hg.type == HYPERGRAPH_TYPE::EDGEWEIGHTS;
	if (edge_weights)
		throw std::runtime_error("Hypergraph in file: " + filename + " has hyperedge weights, which are currently not supported.");
	bool node_weights = hg.type == HYPERGRAPH_TYPE::EDGEANDNODEWEIGHTS || hg.type == HYPERGRAPH_TYPE::NODEWEIGHTS;
	if (node_weights)
		throw std::runtime_error("Hypergraph in file: " + filename + " has node weights, which are currently not supported.");

	nodeid max_net_size = 0;
	for (netid e = 0; e < hg.num_nets; e++) {
		mgetline(f, line);
		std::istringstream iss(line);
		nodeid pin;
		nodeid size_of_net_e = 0;
		while (iss >> pin) {
			if (pin < 1) { throw std::runtime_error("File: " + filename + " has pin id < 1 (in one-based ids)."); }
			else if (pin > hg.num_nodes) { throw std::runtime_error("File: " + filename + " has pin id > number of nodes."); }
			size_of_net_e++;
			hg.pins.push_back(pin-1);
		}
		if (size_of_net_e > hg.num_nodes) { throw std::runtime_error("File: " + filename + " has net with more pins than nodes in the hypergraph."); }
		if (size_of_net_e == 0) { throw std::runtime_error("File: " + filename + " has hyperedge with zero pins."); }


		if (size_of_net_e == 1) { //ignore single pin hyperedges
			hg.pins.pop_back();
		}
		else {
			hg.netsizes.push_back(size_of_net_e);
			max_net_size = std::max(max_net_size, size_of_net_e);
		}
	}
	//std::cout << "max net size " << max_net_size << std::endl;
	//ignore node weights
	f.close();
}

