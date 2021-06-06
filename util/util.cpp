#include "util/util.hpp"

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

bool starts_with(std::string const &value, std::string const &prefix) {
    if (prefix.size() > value.size()) {
        return false;
    }
    return std::equal(prefix.begin(), prefix.end(), value.begin());
}


bool ends_with(std::string const &value, std::string const &ending) {
    if (ending.size() > value.size()) {
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
    for (uint32_t r = 0; r < values.size(); ++r) {
        std::copy(values[r].begin(), values[r].end(), result.row(r));
    }
    return result;
}

uint8_t chr_to_idx(const std::string &chromosome) {
    if (chromosome == "X") {
        return 22;
    } else  if (chromosome == "Y") {
        return 23;
    }
    if (!std::all_of(chromosome.begin(), chromosome.end(), ::isdigit)) {
        return 255;
    }
    return std::stoi(chromosome) - 1;
}

std::vector<std::vector<uint32_t>> read_positions(const std::string &file) {
    if (!std::filesystem::exists(file)) {
        logger()->error("Could not file positions file: {}", file);
        std::exit(1);
    }
    logger()->info("Reading valid positions...");

    std::vector<std::vector<uint32_t>> result;
    std::string line;
    std::ifstream f(file);
    uint32_t pos_count = 0;
    while (std::getline(f, line)) {
        if (starts_with(line, "#")) {
            continue;
        }
        std::istringstream iss(line);

        std::string chromosome_str;
        std::getline(iss, chromosome_str, '\t');
        uint32_t chromosome = chr_to_idx(chromosome_str);
        if (chromosome > 23) {
            continue;
        }
        std::string pos_str;
        std::getline(iss, pos_str, '\t');
        uint32_t pos = std::stoul(pos_str);
        if (chromosome >= result.size()) {
            result.resize(chromosome + 1);
        }
        result[chromosome].push_back(pos);
        pos_count++;
    }
    for (auto& v : result) {
        std::sort(v.begin(), v.end());
    }
    logger()->info("Found a total of  {} positions in {} chromosomes", pos_count, result.size());
    return result;
}

uint32_t chromosome_to_id(const std::string &chromosome) {
    char *p;
    uint32_t converted = strtol(chromosome.c_str(), &p, 10);
    if (*p) {
        if (chromosome != "X" && chromosome != "Y") {
            logger()->error("Invalid chromosome: {}. Must be 1..22, X, Y", chromosome);
            std::exit(1);
        }
        return (chromosome == "X") ? 22 : 23;
    } else {
        return converted - 1;
    }
}


std::string id_to_chromosome(uint32_t chr_id) {
    if (chr_id < 22) {
        return std::to_string(chr_id + 1);
    }
    if (chr_id == 22)
        return "X";
    assert(chr_id == 23);
    return "Y";
}

