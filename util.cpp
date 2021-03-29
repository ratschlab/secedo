#include "util.hpp"

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

std::vector<double> double_split(const std::string &s, char c) {
    std::string segment;
    std::vector<double> result;
    std::stringstream in(s);

    while (std::getline(in, segment, c)) {
        result.push_back(std::stod(segment));
    }
    return result;
}

std::vector<std::filesystem::path> get_files(const std::filesystem::path &path,
                                             const std::string &extension) {
    assert(extension[0] == '.');
    std::vector<std::filesystem::path> result;
    for (auto &p : std::filesystem::recursive_directory_iterator(path)) {
        if (p.path().extension() == extension) {
            result.push_back(p.path());
        }
    }
    return result;
}
