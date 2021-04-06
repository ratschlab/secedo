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

bool ends_with(std::string const &value, std::string const &ending)
{
    if (ending.size() > value.size())
    {
        return false;
    }
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

Matd read_mat(const std::string &name) {
    if (!std::filesystem::exists(name)) {
        logger()->error("File {} does not exist.", name);
        std::exit(1);
    }
    std::ifstream in(name);
    std::string line;
    std::vector<std::vector<double>> values;

    while (std::getline(in, line)) {
        values.push_back(double_split(line, ','));
        if (values.size() > 1 && values.back().size() != values[values.size() - 2].size()) {
            logger()->error("Invalid input file: {}. Line {} has size {}, line {} has size {}",
                            name, values.size() - 2, values[values.size() - 2].size(),
                            values.size() - 1, values.back().size());
            std::exit(1);
        }
    }
    if (values.empty()) {
        return Matd();
    }
    Matd result(values.size(), values[0].size());
    for(uint32_t r = 0; r < values.size(); ++r) {
        std::copy(values[r].begin(), values[r].end(), result.row(r));
    }
    return result;
}

