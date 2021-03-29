#pragma once

#include "expectation_maximization.hpp"
#include "util.hpp"


#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

constexpr uint8_t NO_GENOTYPE = std::numeric_limits<uint8_t>::max();

/**
 * Find the most likely genotype, if we observe bases counts
 * @parm nBases (array of length 4: number of As, Cs, etc.)
 * @param heteroPrior the prior on heterozygous genotype
 * @param theta sequencing error rate
 */
uint8_t
mostLikelyGenotype(const std::array<uint32_t, 4> &nBases, double heteroPrior, double theta) {
    double logTheta = std::log(theta / 3);
    double logOneMinusTheta = std::log(1 - theta);
    double logHalfMinusTheta = std::log(0.5 - theta / 3);

    std::vector<uint32_t> idx = argsort(nBases.begin(), nBases.end());
    // coverage for the given position; // we require coverage at least 9
    uint32_t cov = sum(nBases.begin(), nBases.end());
    if (cov < 9) {
        return NO_GENOTYPE;
    }

    // probability of the most likely homozygous genotype
    double logProb_homo = idx[3] * logOneMinusTheta + (cov - idx[3]) * logTheta;
    // probability of the most likely heterozygous genotype
    double logProb_hetero = (idx[2] + idx[3]) * logHalfMinusTheta + (idx[0] + idx[1]) * logTheta
            + std::log(heteroPrior);

    // if logProb_homo == logProb_hetero, we are not able to decide
    if (logProb_homo == logProb_hetero) {
        return NO_GENOTYPE;
    }

    // homozygous genotype
    if (logProb_homo > logProb_hetero) {
        // check that the genotype is unique
        if (idx[2] == idx[3]) {
            return NO_GENOTYPE;
        }
        char base = idx[3];
        return 4 * base + base;
    }

    // heterozygous genotype
    if (idx[2] == idx[1]) // check uniqueness
        return NO_GENOTYPE;
    return 4 * idx[3] + idx[2];
}

void variant_calling(const std::vector<std::vector<PosData>> &pos_data,
                     const std::vector<double> cell_labels,
                     double hetero_prior,
                     double theta,
                     const std::filesystem::path &out_dir) {
    std::ofstream f(out_dir / "variant");
    uint32_t chromosome = 1;
    for (const std::vector<PosData> &chromosome_data : pos_data) {
        for (const PosData &pd : chromosome_data) {
            std::array<uint32_t, 4> nbases_a;
            std::array<uint32_t, 4> nbases_b;
            for (const CellData &cd : pd.cells_data) {
                if (cell_labels[cd.cell_id] <= 0.05) { // first cluster
                    nbases_a[cd.base]++;
                } else if (cell_labels[cd.cell_id] >= 0.95) {
                    nbases_b[cd.base]++;
                }
            }

            uint8_t genotype_a = mostLikelyGenotype(nbases_a, hetero_prior, theta);
            uint8_t genotype_b = mostLikelyGenotype(nbases_b, hetero_prior, theta);
            // write position to file if different
            if (genotype_a != genotype_b && genotype_a != NO_GENOTYPE
                && genotype_b != NO_GENOTYPE) {
                f << chromosome << pd.position << nbases_a << genotype_a << nbases_b << genotype_b
                  << std::endl;
            }
        }
        chromosome++;
    }
}
