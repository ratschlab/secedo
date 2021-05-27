#include "pileup_reader.hpp"

#include "util/logger.hpp"
#include "util/util.hpp"

#include <filesystem>
#include <iostream>
#include <vector>

constexpr bool simplify = true;

std::tuple<std::vector<PosData>, uint16_t, uint32_t>
read_pileup_text(const std::string fname,
                 const std::vector<uint16_t> &id_to_group,
                 const std::function<void(uint64_t)> &progress,
                 uint32_t max_coverage,
                 const std::vector<uint32_t> &positions) {
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
    std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> id_stats;
    uint32_t id_count = 0;
    uint64_t read_bytes = 0;
    uint32_t pos_idx = 0;
    // process position by position
    while (std::getline(f, line)) {
        read_bytes = (line.size() + 1);
        progress(read_bytes);

        // parse the data, split by tabs
        std::vector<std::string> splitLine = split(line, '\t');

        // 2nd column: position in chromosome
        uint32_t position = std::stoll(splitLine[1]);
        // 4th column: bases/reads for the given position
        std::string bases = splitLine[3];
        // 5th column: cell ids the cell from which the reads are coming, comma separated
        std::vector<uint16_t> cell_ids = int_split<uint16_t>(splitLine[4], ',');

        if (cell_ids.size() > max_coverage) {
            continue;
        }

        if (!positions.empty()) {
            while (pos_idx < positions.size() && positions[pos_idx] < position) {
                pos_idx++;
            }
            if (pos_idx == positions.size()) {
                break; // all valid positions were read
            }
            if (positions[pos_idx] > position) {
                continue;
            }
            assert(positions[pos_idx] == position);
        }

        for (const uint16_t id : cell_ids) {
            all_cell_ids.insert(id);
            all_cell_ids_grouped.insert(id_to_group[id]);
        }

        // 6th column: read ids, comma separated. Each single cell sequence has its own read id
        std::vector<std::string> read_ids = split(splitLine[5], ',');
        std::vector<uint32_t> read_ids_int(read_ids.size());

        for (uint32_t j = 0; j < read_ids.size(); ++j) {
            if (id_map.find(read_ids[j]) == id_map.end()) {
                id_map[read_ids[j]] = id_count++;
            }
            read_ids_int[j] = id_map[read_ids[j]];
        }

        std::vector<uint16_t> cell_ids_and_bases_grouped;
        std::vector<uint16_t> cell_ids_and_bases;
        for (uint32_t j = 0; j < bases.size(); ++j) {
            if (cell_ids[j] > id_to_group.size()) {
                logger()->error(
                        "Cell id {} is too large. Increase --max_cell_count if using the default "
                        "mapping, or fix the mapping in --merge_file",
                        id_to_group.size());
                std::exit(1);
            }

            uint16_t cell_id_base = id_to_group.at(cell_ids[j]) << 2 | CharToInt[(uint8_t)bases[j]];
            cell_ids_and_bases_grouped.push_back(cell_id_base);
            cell_id_base = cell_ids[j] << 2 | CharToInt[(uint8_t)bases[j]];
            cell_ids_and_bases.push_back(cell_id_base);
        }
        result.push_back({ position, read_ids_int, cell_ids_and_bases_grouped });

        for (uint32_t j = 0; j < read_ids.size(); ++j) {
            if (id_stats.find(read_ids_int[j]) == id_stats.end()) {
                id_stats[read_ids_int[j]] = { position, position };
            }
            id_stats[read_ids_int[j]] = { id_stats[read_ids_int[j]].first, position };
        }

        if (simplify) {
            uint16_t coverage = read_ids_int.size();
            out_bin.write(reinterpret_cast<char *>(&position), sizeof(position));
            out_bin.write(reinterpret_cast<char *>(&coverage), sizeof(coverage));
            out_bin.write(reinterpret_cast<char *>(read_ids_int.data()),
                          read_ids_int.size() * sizeof(read_ids_int[0]));
            out_bin.write(reinterpret_cast<char *>(cell_ids_and_bases.data()),
                          cell_ids_and_bases.size() * sizeof(cell_ids_and_bases[0]));
        }
    }

    uint32_t max_length = 0;
    uint32_t max_id = 0;
    for (const auto &[k, v] : id_stats) {
        uint32_t len = v.second - v.first;
        if (max_length < len) {
            max_length = len;
            max_id = k;
        }
    }
    logger()->trace(
            "{}: found {} cell ids, {} after grouping, {} reads. "
            "Longest fragment is {} with {} bases",
            fname, all_cell_ids.size(), all_cell_ids_grouped.size(), id_stats.size(), max_id,
            max_length);

    return { result, all_cell_ids.size(), max_length };
}

std::tuple<std::vector<PosData>, uint16_t, uint32_t>
read_pileup_bin(const std::string fname,
                const std::vector<uint16_t> &id_to_group,
                const std::function<void(uint64_t)> &progress,
                uint32_t max_coverage,
                const std::vector<uint32_t> &positions) {
    std::vector<PosData> result;

    if (!std::filesystem::exists(fname)) {
        spdlog::error("File {} does not exist or is not readable.", fname);
        std::exit(1);
    }

    std::ifstream f(fname, std::ios::binary);
    std::string line;
    uint16_t max_cell_id = 0;
    uint16_t max_cell_id_grouped = 0;

    // maps read ids to a (first,last) pair (in order to find the longest DNA fragment)
    std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> id_stats;

    uint64_t read_bytes = 0;
    uint64_t reported_bytes = 0;
    uint32_t pos_idx = 0;
    uint32_t i = 0;
    while (f.good()) {
        uint32_t position;
        uint16_t coverage;
        f.read(reinterpret_cast<char *>(&position), sizeof(position));
        if (!f.good()) {
            break;
        }
        f.read(reinterpret_cast<char *>(&coverage), sizeof(coverage));

        std::vector<uint32_t> read_ids(coverage);
        std::vector<uint16_t> cell_ids_and_bases(coverage);

        f.read(reinterpret_cast<char *>(read_ids.data()), coverage * sizeof(read_ids[0]));
        f.read(reinterpret_cast<char *>(cell_ids_and_bases.data()),
               coverage * sizeof(cell_ids_and_bases[0]));

        if (i & 4095) {
            // report progress
            read_bytes = f.tellg();
            if ((int64_t)read_bytes == -1) { // end of file
                read_bytes = std::filesystem::file_size(fname);
            }
            progress(read_bytes - reported_bytes);
            reported_bytes = read_bytes;
        }
        i++;
        if (coverage > max_coverage) {
            continue;
        }

        if (!positions.empty()) {
            while (pos_idx < positions.size() && positions[pos_idx] < position) {
                pos_idx++;
            }
            if (pos_idx == positions.size()) {
                break; // all valid positions were read
            }
            if (positions[pos_idx] > position) {
                continue;
            }
            assert(positions[pos_idx] == position);
        }

        for (uint32_t read_idx = 0; read_idx < read_ids.size(); ++read_idx) {
            uint16_t cell_id = cell_ids_and_bases[read_idx] >> 2;
            max_cell_id = std::max(max_cell_id, cell_id);
            if (cell_id >= id_to_group.size()) {
                logger()->error(
                        "Cell id {} is too large. Increase --max_cell_count if using the default "
                        "mapping, or fix the mapping in --merge_file",
                        id_to_group.size());
                std::exit(1);
            }
            uint8_t base = cell_ids_and_bases[read_idx] & 3;
            cell_ids_and_bases[read_idx] = id_to_group[cell_id] << 2 | base;
            max_cell_id_grouped = std::max(max_cell_id_grouped, id_to_group[cell_id]);

            uint32_t read_id = read_ids[read_idx];
            auto it = id_stats.find(read_id);
            if (it != id_stats.end()) {
                it->second = { it->second.first, position };
            } else {
                id_stats[read_id] = { position, position };
            }
        }

        result.emplace_back(position, std::move(read_ids), std::move(cell_ids_and_bases));
    }

    max_cell_id++;
    uint32_t max_length = 0;
    uint32_t max_id = 0;
    for (const auto &[k, v] : id_stats) {
        uint32_t len = v.second - v.first;
        if (max_length < len) {
            max_length = len;
            max_id = k;
        }
    }
    logger()->trace(
            "{}: found {} cell ids, {} after grouping, {} reads. "
            "Longest fragment is {} with {} bases",
            fname, max_cell_id, max_cell_id_grouped, id_stats.size(), max_id, max_length);

    return { result, max_cell_id, max_length };
}

std::tuple<std::vector<PosData>, uint16_t, uint32_t>
read_pileup(const std::string fname,
            const std::vector<uint16_t> &id_to_group,
            const std::function<void(uint64_t)> &progress,
            uint32_t max_coverage,
            const std::vector<uint32_t> positions) {
    return ends_with(fname, ".bin")
            ? read_pileup_bin(fname, id_to_group, progress, max_coverage, positions)
            : read_pileup_text(fname, id_to_group, progress, max_coverage, positions);
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
