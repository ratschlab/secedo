#include "pileup_reader.hpp"

#include "util.hpp"

#include <spdlog/spdlog.h>

#include <filesystem>
#include <vector>

constexpr bool simplify = false;

std::pair<std::vector<PosData>, std::unordered_set<uint32_t>> read_pileup(const std::string fname) {
    std::vector<PosData> result;

    std::ofstream out(fname + ".simplified");

    if (!std::filesystem::exists(fname)) {
        spdlog::error("File {} does not exist or is not readable.", fname);
        std::exit(1);
    }

    std::ifstream f(fname);
    std::string line;
    std::unordered_set<uint32_t> all_cell_ids;
    std::unordered_map<std::string, uint32_t> id_map;
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
        if (simplify) {
            for (uint32_t j = 0; j < read_ids.size(); ++j) {
                if (!id_map.contains(read_ids[j])) {
                    id_map[read_ids[j]] = id_count++;
                }
                read_ids[j] = std::to_string(id_map[read_ids[j]]);
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
        }
    }
    return { result, all_cell_ids };
}
