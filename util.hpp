#pragma once

#include <cstdint>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

constexpr uint8_t CharToInt[128]
        = { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };

template <typename T>
using Mat = std::vector<std::vector<T>>;
using Matd = Mat<double>;
using Mat32u = Mat<uint32_t>;
using Mat64u = Mat<uint64_t>;

template <typename T>
Mat<T> newMat(uint32_t l, uint32_t c, T v = 0) {
    return Mat<T>(l, std::vector<T>(c, v));
}


/** Read an entire file into a string */
std::string read_file(const std::string &fname);

/**
 * Write a vector to a file, with elements separated by comma.
 */
template <typename T>
void write_vec(const std::string &name, const std::vector<T> &vec) {
    std::ofstream out(name);
    for (const auto v : vec) {
        out << v << ", ";
    }
    out << std::endl;
}

/**
 * Write a matrix line by line to a file, with elements separated by comma.
 */
template <typename T>
void write_mat(const std::string &name, const Mat<T> &mat) {
    std::ofstream out(name);
    for (const auto &vec : mat) {
        for (const auto v : vec) {
            out << v << ", ";
        }
        out << std::endl;
    }
}

/**
 * Split a string by character c.
 */
std::vector<std::string> split(const std::string &s, char c);

/**
 * Split a string by character c into integer components.
 */
std::vector<uint32_t> int_split(const std::string &s, char c);

/**
 * Split a string by character c into integer components.
 */
std::vector<double> double_split(const std::string &s, char c);

template <typename T>
T sum(const std::vector<T> &vec) {
    return std::accumulate(vec.begin(), vec.end(), T(0));
}
