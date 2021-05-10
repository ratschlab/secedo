#pragma once

#include <cstdint>
#include <string>
#include <vector>

/**
 * One read for a given cell at a given (unknown) position.
 */
struct CellData {
    uint32_t read_id;
    inline uint16_t cell_id() const {
        return cell_id_base >> 2;
    }
    inline uint8_t base() const {
        return cell_id_base & 3;
    }
    /** First 14 bits are the cell id, last 2 are the base */
    uint16_t cell_id_base;
};

/**
 * All the reads from all the cells for a given position.
 */
struct PosData {
    uint32_t position;
    std::vector<CellData> cells_data;
};
