#pragma once

#include "sequenced_data.hpp"

#include <string>
#include <unordered_set>

std::pair<std::vector<PosData>, std::unordered_set<uint32_t>> read_pileup(const std::string fname);
