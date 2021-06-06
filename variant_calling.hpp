#pragma once

#include "expectation_maximization.hpp"
#include "util/util.hpp"


#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

constexpr uint8_t NO_GENOTYPE = std::numeric_limits<uint8_t>::max();

/**
 * Reads the next chromosome from a FASTA file that contains pairs of (maternal,paternal) contigs
 */
void get_next_chromosome(std::ifstream &fasta_file, std::vector<uint8_t> *chr_data);

/**
 * Calls the most likely variant for each position for each of the clusters given by #clusters.
 * @param pos_data pooled reads for each chromosome at each position; chromosomes are indexed from
 * 0 to 23, with 22 representing chr X and 23 chr Y.
 * @param clusters the cluster to which each cell belongs to
 * @param reference_genome the genome to call the variants against
 * @param hetero_prior the probability that a locus is heterozygous
 * @param theta sequencing error rate
 * @param out_dir location where to write the VCF files that describe the variants
 */
void variant_calling(const std::vector<std::vector<PosData>> &pos_data,
                     const std::vector<uint16_t> &clusters,
                     const std::string &reference_genome,
                     double hetero_prior,
                     double theta,
                     const std::filesystem::path &out_dir);
