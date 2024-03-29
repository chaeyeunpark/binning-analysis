#include <charconv>
#include <iostream>
#include <random>
#include <span>
#include <vector>

#include "LogBinning.hpp"

double energy(std::span<const int32_t> spins) {
	double e = 0;
	for(size_t i = 0; i < spins.size() - 1; i++) {
		e += -spins[i] * spins[i+1];
	}
	return e;
}

double energy(size_t nspins, uint64_t spins) {
	double e = 0;
	for(size_t i = 0; i < nspins - 1; i++) {
		int32_t spin_i = 1-2*static_cast<int32_t>((spins >> i) & 1);
		int32_t spin_j = 1-2*static_cast<int32_t>((spins >> (i+1)) & 1);
		e += - spin_i*spin_j;
	}
	return e;
}

double ediff(std::span<const int32_t> spins, size_t flip_idx) {
	if (flip_idx == 0) {
		return 2.0*spins[flip_idx]*spins[flip_idx+1];
	}
	else if (flip_idx == spins.size() - 1) {
		return 2.0*spins[flip_idx-1]*spins[flip_idx];
	}
	else {
		return 2.0*spins[flip_idx-1]*spins[flip_idx] + 2.0*spins[flip_idx]*spins[flip_idx+1];
	}
}

int main(int argc, char* argv[]) {
	if(argc != 2) {
		std::cout << "Usage: " << argv[0] << " [N]\n";
		return 1;
	}

	const size_t nspins = [&]() {
        std::string_view s = argv[1];
        size_t res{};
        std::from_chars(std::begin(s), std::end(s), res);
        return res;
    }();

	std::mt19937_64 re{1337};
	std::uniform_int_distribution udist(0, 1);
	std::uniform_real_distribution<double> pdist(0.0, 1.0);
	std::uniform_int_distribution<size_t> flip_dist(0, nspins - 1);

	for(const double beta: {0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0}) {
	// for(const double beta: {10.0}) {
		std::vector<int32_t> spins;
		spins.reserve(nspins);

		for(size_t i = 0; i < nspins; i++) {
			spins.push_back(1-2*udist(re));
		}

		LogBinning<double> log_binning(14);
		for(size_t i = 0; i < (1U << 20); i++) {
			double e = energy(spins);
			log_binning.append_data(e);

			for(size_t k = 0; k < 2; k++) {
				size_t flip_idx = flip_dist(re);
				double ed = ediff(spins, flip_idx);
				if(pdist(re) < std::exp(-beta*ed)) {
					spins[flip_idx] *= -1;
				}
			}
		}

		auto mean = log_binning.mean();
		auto corrs = log_binning.vars();

		std::cout << beta << "\t" << mean << "\t";
		for(auto x: corrs) {
			std::cout << x << "\t";
		}
		std::cout << std::endl;

		double Z = 0.0;
		double emean = 0.0;
		double esqr = 0.0;
		for(size_t i = 0; i < (static_cast<size_t>(1U) << nspins); i++) {
			double e = energy(nspins, i);
			double partition = std::exp(-beta*e);
			emean += e * partition;
			esqr += e * e * partition;
			Z += partition;
		}
		emean /= Z;
		esqr /= Z;

		// For uncorrelated sample
		std::cout << emean << "\t" << (esqr - emean*emean) << std::endl;
	}

	return 0;
}
