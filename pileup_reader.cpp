#include "pileup_reader.hpp"

#include "logger.hpp"
#include "util.hpp"

#include <filesystem>
#include <iostream>
#include <vector>

constexpr bool simplify = true;

std::tuple<std::vector<PosData>, std::unordered_set<uint32_t>, uint32_t>
read_pileup_text(const std::string fname) {
    std::vector<PosData> result;

    std::ofstream out(fname + ".simplified");
    std::ofstream out_bin(fname + ".bin", std::ios::binary);

    if (!std::filesystem::exists(fname)) {
        spdlog::error("File {} does not exist or is not readable.", fname);
        std::exit(1);
    }

    std::ifstream f(fname);
    std::string line;
    std::unordered_set<uint32_t> all_cell_ids;
    std::unordered_map<std::string, uint32_t> id_map;
    std::unordered_map<std::string, std::pair<uint32_t, uint32_t>> id_stats;
    uint32_t id_count = 0;
    // process position by position
    while (std::getline(f, line)) {
        // parse the data, split by tabs
        std::vector<std::string> splitLine = split(line, '\t');

        // 2nd column: position in chromosome
        uint64_t position = std::stoll(splitLine[1]);
        // 4th column: bases/reads for the given position
        std::string bases = splitLine[3];
        // 5th column: cell ids the cell from which the reads are coming, comma separated
        std::vector<uint16_t> cell_ids = int_split<uint16_t>(splitLine[4], ',');
        std::copy(cell_ids.begin(), cell_ids.end(),
                  std::inserter(all_cell_ids, all_cell_ids.end()));

        // 6th column: read ids, comma separated. Each single cell sequence has its own read id
        std::vector<std::string> read_ids = split(splitLine[5], ',');
        std::vector<uint32_t> read_ids_int(read_ids.size());
        if (simplify) {
            for (uint32_t j = 0; j < read_ids.size(); ++j) {
                if (!id_map.contains(read_ids[j])) {
                    id_map[read_ids[j]] = id_count++;
                    id_stats[read_ids[j]] = { position, position };
                }
                id_stats[read_ids[j]] = { id_stats[read_ids[j]].first, position };
                read_ids_int[j] = id_map[read_ids[j]];
            }
        }

        std::vector<CellData> cells_data;
        for (uint32_t j = 0; j < bases.size(); ++j) {
            cells_data.push_back({ read_ids[j], cell_ids[j], CharToInt[(uint8_t)bases[j]] });
        }
        result.push_back({ position, cells_data });

        if (simplify) {
            for (uint32_t i = 0; i < splitLine.size(); ++i) {
                if (i == 5) {
                    std::string str_ids = to_string(read_ids);
                    out << str_ids.substr(1, str_ids.size() - 2);
                } else {
                    out << splitLine[i];
                }
                out << '\t';
            }
            out << std::endl;

            uint8_t chromosome = std::stoi(splitLine[0]);
            out_bin.write(reinterpret_cast<char *>(&chromosome), 1);
            out_bin.write(reinterpret_cast<char *>(&position), sizeof(position));
            uint16_t coverage = bases.size();
            out_bin.write(reinterpret_cast<char *>(&coverage), sizeof(coverage));
            out_bin.write(reinterpret_cast<char *>(bases.data()), bases.size());
            out_bin.write(reinterpret_cast<char *>(cell_ids.data()),
                          cell_ids.size() * sizeof(cell_ids.data()[0]));
            out_bin.write(reinterpret_cast<char *>(read_ids_int.data()),
                          read_ids_int.size() * sizeof(read_ids_int.data()[0]));
        }
    }

    uint32_t max_length = 0;
    for (const auto &[k, v] : id_stats) {
        max_length = std::max(max_length, v.second - v.first);
    }
    logger()->trace("{}: found {} cell ids, {} reads. Longest read has {} bases", fname,
                    all_cell_ids.size(), id_stats.size(), max_length);

    return { result, all_cell_ids, max_length };
}

std::tuple<std::vector<PosData>, std::unordered_set<uint32_t>, uint32_t>
read_pileup_bin(const std::string fname) {
    std::vector<PosData> result;

    if (!std::filesystem::exists(fname)) {
        spdlog::error("File {} does not exist or is not readable.", fname);
        std::exit(1);
    }

    std::ifstream f(fname, std::ios::binary);
    std::string line;
    std::unordered_set<uint32_t> all_cell_ids;

    std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> id_stats;

    while (f.good()) {
        uint8_t chromosome;
        uint64_t position;
        uint16_t coverage;
        f.read(reinterpret_cast<char *>(&chromosome), 1);
        if (!f.good()) {
            break;
        }
        f.read(reinterpret_cast<char *>(&position), sizeof(position));
        f.read(reinterpret_cast<char *>(&coverage), sizeof(coverage));

        std::vector<char> bases(coverage);
        std::vector<uint16_t> cell_ids(coverage);
        std::vector<uint32_t> read_ids(coverage);

        f.read(bases.data(), bases.size());
        f.read(reinterpret_cast<char *>(cell_ids.data()), cell_ids.size() * sizeof(cell_ids[0]));
        f.read(reinterpret_cast<char *>(read_ids.data()), read_ids.size() * sizeof(read_ids[0]));

        std::copy(cell_ids.begin(), cell_ids.end(),
                  std::inserter(all_cell_ids, all_cell_ids.end()));

        std::vector<CellData> cells_data;
        for (uint32_t j = 0; j < bases.size(); ++j) {
            if (id_stats.contains(read_ids[j])) {
                id_stats[read_ids[j]] = { id_stats[read_ids[j]].first, position };
            } else {
                id_stats[read_ids[j]] = { position, position };
            }
            cells_data.push_back(
                    { std::to_string(read_ids[j]), cell_ids[j], CharToInt[(uint8_t)bases[j]] });
        }
        result.push_back({ position, cells_data });
    }

    uint32_t max_length = 0;
    for (const auto &[k, v] : id_stats) {
        max_length = std::max(max_length, v.second - v.first);
    }

    logger()->trace("{}: found {} cell ids, {} reads. Longest read has {} bases", fname,
                    all_cell_ids.size(), id_stats.size(), max_length);

    return { result, all_cell_ids, max_length };
}

std::tuple<std::vector<PosData>, std::unordered_set<uint32_t>, uint32_t>
read_pileup(const std::string fname) {
    return fname.ends_with(".bin") ? read_pileup_bin(fname) : read_pileup_text(fname);
}
