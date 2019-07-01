#include <stdexcept>
#include "extern/patoh_interface.h"
#include "io/hmetisreader.h"

int main(int argc, const char* argv[]) {
	if (argc != 5) {
		throw std::runtime_error("Must provide four arguments. Hypergraph, epsilon, seed, PaToH-preset (Q(uality) or D(default)).");
	}

	std::string path_hg = argv[1];
	std::string str_eps = argv[2];
	std::string str_seed = argv[3];
	std::string preset = argv[4];
	double epsilon = std::stod(str_eps);
	int seed = std::stoi(str_seed);
	std::vector<int> seeds = {seed};

	hyper::Hypergraph hg = hyper::HMetisReader::read(path_hg);
	auto t_begin = hyper::time_now();
	std::vector<hyper::Partition> partitions = PaToHInterface::bisectWithPatoh(hg, seeds, epsilon, preset);
	auto t_end = hyper::time_now();

	nodeid perf_balance = hg.totalNodeWeight() / 2;
	auto desired_smb = hyper::Metrics::smallerBlockSize(hg.totalNodeWeight(), epsilon);
	if (epsilon == 0.0 || desired_smb > perf_balance) desired_smb = perf_balance;	//just being cautios. in case of rounding errors

	for (auto& p : partitions) {
		nodeid a0 = p.getSetSize(0);
		nodeid b0 = p.getSetSize(1);
		auto[achieved_smb, ma] = std::minmax(a0, b0);
		double achieved_eps = hyper::Metrics::imbalance(hg.totalNodeWeight(), achieved_smb);

		hyper::flow_t cut = p.computeCut(hg);
		//Format: algoname,hg,eps,achieved_eps,seed,cut,runtime
		std::cout
					<< "PaToH-" << preset << ","
					<< path_hg.substr(path_hg.find_last_of('/') + 1) << ","
					<< epsilon << ","
					<< achieved_eps << ","
					<< str_seed << ","
					<< cut << ","
					<< hyper::inSeconds(hyper::duration(t_end - t_begin)).count()
					<< '\n'
		;
	}
}