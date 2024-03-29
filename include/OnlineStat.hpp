#pragma once
#include <cstddef>

template<class T> class OnlineStat
{
    // Number of samples in this bin
    size_t n_data_;

    // current mean
    T mean_;

    // Running sum of squares of differences from the current mean
    T m2_;

public:
    explicit OnlineStat() : n_data_{0}, mean_{}, m2_{} {
	}

    // Adding data to this bin
    void append_data(const T& data)
    {
        ++n_data_;
        const T delta = data - mean_;
        mean_ += delta / double(n_data_);

        const T delta2 = data - mean_;
        m2_ += delta * delta2;
    }

    void merge_with(const OnlineStat<T>& other)
    {
        n_data_ += other.n_data_;
        const T delta = other.mean_ - mean_;
        mean_ += delta * other.n_data / double(n_data_);
        m2_ += other.m2_ + delta * delta * other.n_data_ * (1 - other.n_data / double(n_data_));
    }

    size_t n_data() const { return n_data_; }

    T mean() const { return mean_; }

    T var() const { return m2_ / double(n_data_ - 1); }

    T error_of_mean() const { return sqrt(m2_ / double(n_data_ * (n_data_ - 1))); }

    void reset()
    {
        n_data_ = 0;
        mean_ = T{};
        m2_ = T{};
    }
};
