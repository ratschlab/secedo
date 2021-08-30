#pragma once

#include <cstdint>
#include <filesystem>
#include <vector>

#include "sequenced_data.hpp"

/**
 * Creates pileup files for a given chromosome from the input BAM files.
 * @param bam_files a list of BAM files that are going to be piled up
 * (https://en.wikipedia.org/wiki/Pileup_format)
 * @param out_pileup path where the resulting pileup file will be written
 * @param write_text_file if true, a text pileup file will be generated for debugging purposes;
 * otherwise the pileup information will be written in our own binary format to save disk space
 * @param chromosome_id the chromosome (1-23) for which information will be piled up
 * @param max_coverage remove positions that have coverage > max_coverage as noise due to
 * amplificaiton artifacts (greatly improves results in real datasets)
 * @param min_base_quality only consider positions with the sequencing quality (phred score) higher
 * than the given threshold (thresold 20=1% errors, threshold 30=0.1% errors)
 * @param min_map_quality only consider sequences mapped with at least the given confidence (if
 * min_map_quality<10, it meanse there is a 1 in 10 chance that the string was misaligned)
 *
 * @param num_threads threads to use for processing the information
 * @param min_different minimum number of bases that must be different from the majority for a locus
 * to be included in the pileup
 * @return the pileup information, one entry for each allele
 */
std::vector<PosData> pileup_bams(const std::vector<std::filesystem::path> &bam_files,
                                 const std::filesystem::path &out_pileup,
                                 bool write_text_file,
                                 uint32_t chromosome_id,
                                 uint32_t max_coverage,
                                 uint32_t min_base_quality,
                                 uint32_t min_map_quality,
                                 uint32_t min_alignment_score,
                                 uint32_t num_threads,
                                 uint16_t min_different);
