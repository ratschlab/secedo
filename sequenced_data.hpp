#pragma once

#include <cassert>
#include <cstdint>
#include <string>
#include <vector>

/**
 * All the reads from all the cells for a given position.
 */
struct PosData {
    PosData(uint32_t position,
            std::vector<uint32_t> &&read_ids,
            std::vector<uint16_t> &&cell_ids_bases)
        : position(position),
          read_ids(std::move(read_ids)),
          cell_ids_bases(std::move(cell_ids_bases)) {}

    PosData(uint32_t position,
            const std::vector<uint32_t> &read_ids,
            const std::vector<uint16_t> &cell_ids_bases)
            : position(position),
              read_ids(std::move(read_ids)),
              cell_ids_bases(std::move(cell_ids_bases)) {}

    uint32_t position;
    /** ids of all reads at this position */
    std::vector<uint32_t> read_ids;
    /**
     * Cell ids and bases of all reads at this position. First 14 bits are the cell id, last 2 are
     * the base
     */
    std::vector<uint16_t> cell_ids_bases;

    inline uint16_t cell_id(uint32_t pos) const { return cell_ids_bases[pos] >> 2; }
    inline uint8_t base(uint32_t pos) const { return cell_ids_bases[pos] & 3; }
    inline uint32_t size() const {
        assert(read_ids.size() == cell_ids_bases.size());
        return read_ids.size();
    }

    bool operator==(const PosData &other) const {
        return position == other.position && read_ids == other.read_ids
                && cell_ids_bases == other.cell_ids_bases;
    }
};
