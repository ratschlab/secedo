#pragma once

#include <cstdint>
#include <string>
#include <vector>


/**
 * One read for a given cell at a given (unknown) position.
 */
struct CellData {
    uint32_t read_id;
    uint16_t cell_id;
    uint8_t base;
};

/**
 * All the reads from all the cells for a given position.
 * Note that the position itself is not stored, as its never being used.
 */
struct PosData {
    uint64_t position;
    std::vector<CellData> cells_data;
};
