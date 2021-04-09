#include "pileup_reader.hpp"

#include "logger.hpp"
#include "util.hpp"

#include <filesystem>
#include <iostream>
#include <vector>

constexpr bool simplify = true;

std::tuple<std::vector<PosData>, std::unordered_set<uint32_t>, uint32_t>
read_pileup_text(const std::string fname,
                 const std::vector<uint16_t> &id_to_group,
                 const std::function<void(uint64_t)> &progress) {
    std::vector<PosData> result;

    std::ofstream out_bin(fname + ".bin", std::ios::binary);

    if (!std::filesystem::exists(fname)) {
        spdlog::error("File {} does not exist or is not readable.", fname);
        std::exit(1);
    }

    std::ifstream f(fname);
    std::string line;
    std::unordered_set<uint32_t> all_cell_ids;
    std::unordered_set<uint32_t> all_cell_ids_grouped;

    std::unordered_map<std::string, uint32_t> id_map;
    std::unordered_map<std::string, std::pair<uint32_t, uint32_t>> id_stats;
    uint32_t id_count = 0;
    uint64_t read_bytes = 0;
    // process position by position
    while (std::getline(f, line)) {
        read_bytes = (line.size() + 1);
        progress(read_bytes);

        // parse the data, split by tabs
        std::vector<std::string> splitLine = split(line, '\t');

        // 2nd column: position in chromosome
        uint64_t position = std::stoll(splitLine[1]);
        // 4th column: bases/reads for the given position
        std::string bases = splitLine[3];
        // 5th column: cell ids the cell from which the reads are coming, comma separated
        std::vector<uint16_t> cell_ids = int_split<uint16_t>(splitLine[4], ',');

        for (const uint16_t id : cell_ids) {
            all_cell_ids.insert(id);
            all_cell_ids_grouped.insert(id_to_group[id]);
        }

        // 6th column: read ids, comma separated. Each single cell sequence has its own read id
        std::vector<std::string> read_ids = split(splitLine[5], ',');
        std::vector<uint32_t> read_ids_int(read_ids.size());
        if (simplify) {
            for (uint32_t j = 0; j < read_ids.size(); ++j) {
                if (id_map.find(read_ids[j]) == id_map.end()) {
                    id_map[read_ids[j]] = id_count++;
                    id_stats[read_ids[j]] = { position, position };
                }
                id_stats[read_ids[j]] = { id_stats[read_ids[j]].first, position };
                read_ids_int[j] = id_map[read_ids[j]];
            }
        }

        std::vector<CellData> cells_data;
        for (uint32_t j = 0; j < bases.size(); ++j) {
            if (cell_ids[j] > id_to_group.size()) {
                logger()->error(
                        "Cell id {} is too large. Increase --max_cell_count if using the default "
                        "mapping, or fix the mapping in --merge_file",
                        id_to_group.size());
                std::exit(1);
            }

            cells_data.push_back(
                    { read_ids[j], id_to_group.at(cell_ids[j]), CharToInt[(uint8_t)bases[j]] });
        }
        result.push_back({ position, cells_data });

        if (simplify) {
            uint8_t chromosome = splitLine[0] == "X" ? 23 : std::stoi(splitLine[0]);
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
    logger()->trace("{}: found {} cell ids, {} after grouping, {} reads. Longest read has {} bases",
                    fname, all_cell_ids.size(), all_cell_ids_grouped.size(), id_stats.size(),
                    max_length);

    return { result, all_cell_ids, max_length };
}

std::tuple<std::vector<PosData>, std::unordered_set<uint32_t>, uint32_t>
read_pileup_bin(const std::string fname,
                const std::vector<uint16_t> &id_to_group,
                const std::function<void(uint64_t)> &progress) {
    std::vector<PosData> result;

    if (!std::filesystem::exists(fname)) {
        spdlog::error("File {} does not exist or is not readable.", fname);
        std::exit(1);
    }

    std::ifstream f(fname, std::ios::binary);
    std::string line;
    std::unordered_set<uint32_t> all_cell_ids;
    std::unordered_set<uint32_t> all_cell_ids_grouped;

    std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> id_stats;

    uint64_t read_bytes = 0;
    uint64_t reported_bytes = 0;

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

        // report progress
        read_bytes = f.tellg();
        if ((int64_t )read_bytes == -1 ) { // end of file
            read_bytes = std::filesystem::file_size(fname);
        }
        progress(read_bytes - reported_bytes);
        reported_bytes = read_bytes;

        for (const uint16_t id : cell_ids) {
            all_cell_ids.insert(id);
            all_cell_ids_grouped.insert(id_to_group[id]);
        }

        std::vector<CellData> cells_data;
        for (uint32_t j = 0; j < bases.size(); ++j) {
            if (cell_ids[j] > id_to_group.size()) {
                logger()->error(
                        "Cell id {} is too large. Increase --max_cell_count if using the default "
                        "mapping, or fix the mapping in --merge_file",
                        id_to_group.size());
                std::exit(1);
            }

            if (id_stats.find(read_ids[j]) != id_stats.end()) {
                id_stats[read_ids[j]] = { id_stats[read_ids[j]].first, position };
            } else {
                id_stats[read_ids[j]] = { position, position };
            }
            cells_data.push_back({ std::to_string(read_ids[j]), id_to_group.at(cell_ids[j]),
                                   CharToInt[(uint8_t)bases[j]] });
        }
        result.push_back({ position, cells_data });
    }

    uint32_t max_length = 0;
    for (const auto &[k, v] : id_stats) {
        max_length = std::max(max_length, v.second - v.first);
    }

    logger()->trace("{}: found {} cell ids, {} after grouping, {} reads. Longest read has {} bases",
                    fname, all_cell_ids.size(), all_cell_ids_grouped.size(), id_stats.size(),
                    max_length);

    return { result, all_cell_ids, max_length };
}

std::tuple<std::vector<PosData>, std::unordered_set<uint32_t>, uint32_t>
read_pileup(const std::string fname,
            const std::vector<uint16_t> &id_to_group,
            const std::function<void(uint64_t)> &progress) {
    return ends_with(fname, ".bin") ? read_pileup_bin(fname, id_to_group, progress)
                                    : read_pileup_text(fname, id_to_group, progress);
}


std::vector<uint16_t>
get_grouping(uint16_t merge_count, const std::string &merge_file, uint16_t max_cell_count) {
    assert(merge_count == 1 || merge_file.empty());

    std::vector<uint16_t> cellid_to_group;
    if (!merge_file.empty()) {
        if (!std::filesystem::exists(merge_file)) {
            logger()->error("Cannot find merge file: {}", merge_file);
            std::exit(1);
        }
        cellid_to_group = int_split<uint16_t>(read_file(merge_file), ',');
    } else {
        for (uint16_t i = 0; i < max_cell_count; ++i) {
            cellid_to_group.push_back(i / merge_count);
        }
    }

    return cellid_to_group;
}
