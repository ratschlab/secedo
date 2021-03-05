#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

template <typename T>
using Mat = std::vector<std::vector<T>>;
using Matd = Mat<double>;
using Mat32u = Mat<uint32_t>;

template <typename T>
Mat<T> newMat(uint32_t l, uint32_t c, T v = 0) {
    return Mat<T>(l, std::vector<T>(c, v));
}


/** Read an entire file into a string */
std::string read_file(const std::string &fname);

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
