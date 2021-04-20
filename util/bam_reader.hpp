#pragma once

#include "heap.hpp"
#include "util/logger.hpp"
#include "preprocess.hpp"
#include "sequenced_data.hpp"
#include "util/util.hpp"

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
#include <progress_bar/progress_bar.hpp>


#include <array>
#include <atomic>
#include <filesystem>
#include <string>
#include <vector>

std::vector<PosData> read_bam(const std::vector<std::filesystem::path> &input_files,
                              const std::filesystem::path &outfile,
                              uint32_t chromosome_id,
                              uint32_t max_coverage,
                              uint32_t num_threads,
                              double sequencing_error_rate);
