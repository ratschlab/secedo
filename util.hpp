#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>


/** Read an entire file into a string */
std::string read_file(const std::string &fname);

/**
 * Write a matrix line by line to a file, with elements separated by comma.
 */
template <typename T>
void write_mat(const std::string &name, const std::vector<std::vector<T>> &mat) {
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
