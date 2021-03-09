#include "util.hpp"

#include <sstream>

std::string read_file(const std::string &fname) {
    std::ifstream f(fname);
    std::string str;

    f.seekg(0, std::ios::end);
    str.reserve(f.tellg());
    f.seekg(0, std::ios::beg);

    str.assign((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    return str;
}

std::vector<std::string> split(const std::string &s, char c) {
    std::string segment;
    std::vector<std::string> result;
    std::stringstream in(s);

    while (std::getline(in, segment, c)) {
        result.push_back(segment);
    }
    return result;
}

std::vector<uint32_t> int_split(const std::string &s, char c) {
    std::string segment;
    std::vector<uint32_t> result;
    std::stringstream in(s);

    while (std::getline(in, segment, c)) {
        result.push_back(std::stoi(segment));
    }
    return result;
}
