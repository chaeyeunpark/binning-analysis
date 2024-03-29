#pragma once

#include "OnlineStat.hpp"

#include <vector>

template <typename T>
class LogBinning {
private:
	size_t max_level_;
	std::vector<T> buffer_;
	std::vector<OnlineStat<T>> bins_;

public:
	explicit LogBinning(size_t max_level) : max_level_{max_level}, bins_(max_level+1) {
	}

	void append_data(T data) {
		buffer_.emplace_back(std::move(data));

		if(buffer_.size() == static_cast<size_t>(1U << max_level_)) {
			std::vector<T> buf_curr_level = buffer_;
			for(size_t curr_level = 0; curr_level < max_level_ + 1; curr_level++)
			{
				for(size_t i = 0; i < buf_curr_level.size(); i++) {
					bins_[curr_level].append_data(buf_curr_level[i]);
				}

				std::vector<T> buf_next_level;
				for(size_t i = 0; i < buf_curr_level.size() / 2; i++) {
					buf_next_level.push_back((buf_curr_level[2*i]+buf_curr_level[2*i+1])/2.0);
				}

				buf_curr_level.swap(buf_next_level);
			}

			std::vector<T>{}.swap(buffer_);
		}
	}

	T mean() const {
		return bins_[0].mean();
	}

	std::vector<T> vars() const {
		std::vector<T> res;
		for(size_t i = 0; i < max_level_+1; i++) {
			res.push_back(bins_[i].var());
		}
		return res;
	}

	std::vector<T> error_of_means() const {
		std::vector<T> res;
		for(size_t i = 0; i < max_level_+1; i++) {
			res.push_back(bins_[i].error_of_mean());
		}
		return res;
	}

	std::vector<T> corrs() const {
		std::vector<T> res;
		for(size_t i = 0; i < max_level_ + 1; i++) {
			T corr = 0.5*((bins_[i].var() / bins_[i].n_data())/(bins_[0].var() / bins_[0].n_data()) - 1);
			res.emplace_back(std::move(corr));
		}
		return res;
	}
};
