#pragma once

// inlude only here, so that all files use the same default logger
#include <spdlog/spdlog.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "util/logger.hpp"
#include "util/mat.hpp"

constexpr uint8_t CharToInt[128]
        = { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 5, 1, 5, 5, 5, 2,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };

constexpr char IntToChar[6] = { 'A', 'C', 'G', 'T', 'N', 'N' };

const std::string int_to_chr[24]
        = { "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12",
            "13", "14", "15", "16", "17", "18", "19", "20", "21", "X",  "Y" };

template <typename T>
using Vec2 = std::vector<std::vector<T>>;
using Vec2d = Vec2<double>;
using Vec2u32 = Vec2<uint32_t>;
using Vec2u64 = Vec2<uint64_t>;

template <typename T>
Mat<T> newVec2(uint32_t l, uint32_t c, T v = 0) {
    return Vec2<T>(l, std::vector<T>(c, v));
}


/** Read an entire file into a string */
std::string read_file(const std::string &fname);

/**
 * Write a vector to a file, with elements separated by comma.
 */
template <typename T>
void write_vec(const std::string &name, const std::vector<T> &vec) {
    std::ofstream out(name);
    if (vec.empty()) {
        return;
    }

    for (uint32_t i = 0; i < vec.size() - 1; ++i) {
        out << vec[i] << ",";
    }
    out << vec.back();
    out << std::endl;
}

/**
 * Join a vector's element into comma-separated string.
 */
template <typename T>
std::string join_vec(const std::vector<T> &vec) {
    std::stringstream out;
    if (vec.empty()) {
        return "";
    }

    for (uint32_t i = 0; i < vec.size() - 1; ++i) {
        out << vec[i] << ",";
    }
    out << vec.back();
    out << std::endl;
    return out.str();
}


/**
 * Write a vector to a string, with elements separated by comma.
 */
template <typename T>
std::string to_string(const std::vector<T> &vec) {
    if (vec.empty()) {
        return "[]";
    }
    std::stringstream out;
    out << '[';
    for (uint32_t i = 0; i < vec.size() - 1; ++i) {
        out << vec[i] << ",";
    }
    out << vec.back() << ']';
    return out.str();
}

/**
 * Split a string by character c.
 */
std::vector<std::string> split(const std::string &s, char c);

/**
 * Split a string by character c into integer components.
 */
template <typename T>
std::vector<T> int_split(const std::string &s, char c) {
    std::string segment;
    std::vector<T> result;
    std::stringstream in(s);
    try {
        while (std::getline(in, segment, c)) {
            result.push_back(std::stoll(segment));
        }
    } catch (std::invalid_argument &) {
        logger()->error("Invalid string to split: {}", s);
        std::exit(1);
    }
    return result;
}

/**
 * Split a string by character c into integer components.
 */
std::vector<double> double_split(const std::string &s, char c);

/**
 * Write a matrix line by line to a file, with elements separated by comma.
 */
template <typename T>
void write_mat(const std::string &name, const Mat<T> &mat) {
    std::ofstream out(name);
    if (mat.empty()) {
        return;
    }
    for (uint32_t r = 0; r < mat.rows(); ++r) {
        for (uint32_t i = 0; i < mat.cols() - 1; ++i) {
            out << mat(r, i) << ",";
        }
        out << mat(r, mat.cols() - 1) << std::endl;
    }
}

/**
 * Read a double matrix written with #write_mat.
 */
Matd read_mat(const std::string &name);

template <typename T>
T sum(const std::vector<T> &vec) {
    return std::accumulate(vec.begin(), vec.end(), T(0));
}

template <typename iterator>
typename std::iterator_traits<iterator>::value_type sum(iterator b, iterator e) {
    typename std::iterator_traits<iterator>::value_type z(0);
    return std::accumulate(b, e, z);
}

/**
 * Returns a permutation of {0....a.size()-1}, such that the elements of a are sorted.
 * @param b beginning iterator for data to sort
 * @param e end iterator for data to sort
 * @return permutation of the elements such that they are sorted in ascending order
 */
template <typename iterator>
std::vector<uint32_t> argsort(iterator b, iterator e) {
    std::vector<uint32_t> result(std::distance(b, e));
    std::iota(result.begin(), result.end(), 0);
    std::sort(result.begin(), result.end(), [&b](int left, int right) -> bool {
        // sort indices according to the corresponding element
        return *(b + left) < *(b + right);
    });

    return result;
}

template <typename T, uint32_t C>
std::array<uint32_t, C> argsort(const std::array<T, C> &array) {
    std::array<uint32_t, C> result;
    std::iota(result.begin(), result.end(), 0);
    std::sort(result.begin(), result.end(),
              [&array](int left, int right) -> bool { return array[left] < array[right]; });

    return result;
}

template <class T>
inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
    os << "[";
    for (auto &x : v) {
        os << " " << x;
    }
    os << " ]";
    return os;
}

template <class T, size_t N>
inline std::ostream &operator<<(std::ostream &os, const std::array<T, N> &v) {
    os << "[";
    for (auto &x : v) {
        os << " " << x;
    }
    os << " ]";
    return os;
}

/**
 * Find files with the given extension in the given path (recursive)
 * @param path where to look for files
 * @param extension the extension to look for (must start with '.')
 * @return the list of files in #path that match #extension
 */
std::vector<std::filesystem::path> get_files(const std::filesystem::path &path,
                                             const std::string &extension);

bool ends_with(std::string const &value, std::string const &ending);


std::vector<std::vector<uint32_t>> read_positions(const std::string &file);

/**
 * Converts a human chromosome names (1..22, X, Y) to an id between 0 and 23
 * @param chromosome a one character string representing the chromosome name
 * @return an id between 0 and 23, with 22=X, 23=Y
 */
uint8_t chromosome_to_id(const std::string &chromosome);

/** Converts a chromosome id between 0 and 22 to the chromosome name 1..22, X, Y */
std::string id_to_chromosome(uint32_t chr_id);
