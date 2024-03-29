#include "is_significant.hpp"

#include "util/logger.hpp"

#include <cfenv>
#include <numeric>

// these values were computed using scripts/K.r and correspond to the optimal threshold (acc to
// Youden's J) for a 10-90, 20-80, 30-70, 40-60 and 50-50 split respectively between cells
// Values computed for p_mut=1e-6, p_het=5e-4
std::vector<std::vector<double>> Filter::Ks
        = { { -1.64504967001201, -1.38868450353301, -1.38780664765677, -1.38779600211955,
              -1.3877952855556,  -1.38779524274215, -1.38779524274142, -1.38779524274141,
              -1.3877952427414,  -1.38779524274139, -1.38779524274138, -1.38779524274138,
              -1.38779524274138, -1.38779524274138, -1.38779524274139, -1.38780870444455,
              -1.38780870444455, -1.38780870444455, -1.38780870444455, -1.38780870444455 },
            { -1.56013904495168, -1.38819451352203, -1.38781438946096, -1.38779659244035,
              -1.38779537799054, -1.3877952484612,  -1.3877952427842,  -1.38779524274906,
              -1.3877952427457,  -1.38779524274275, -1.38780870444458, -1.38780870444459,
              -1.38780870444459, -1.38780870444469, -1.38780870444459, -1.42736056742577,
              -1.42736056742575, -6.19144172018466, -6.19144172018466, -14.1885779508362 },
            { -1.47780038365618, -1.3885722463397,  -1.38781428162649, -1.3877984410546,
              -1.38779548312685, -1.3877952855556,  -1.38779524455204, -1.38779524331456,
              -1.38780873675669, -1.38780870687333, -1.42737804009806, -6.19144172131432,
              -14.1885779508648, -6.1914418045659,  -30.1993093269287, -30.1993093268559,
              -30.1993093268539, -54.2154105288032, -62.2207775961199, -46.2100434614866 },
            { -1.47780038365618, -1.38868450353301, -1.38782829051844, -1.3877984410546,
              -1.38779625512927, -1.38779556321717, -1.38780972304588, -1.3878087226245,
              -6.21747860711653, -22.1939432034943, -14.1886670526002, -22.1939422903721,
              -46.2100434614866, -54.2154105288069, -70.2261446634366, -62.2207775961199,
              -86.2368787980699, -110.25298000002,  -118.258347067337, -102.247612932703 },
            { -1.52859626647315, -1.38967447346712, -1.38787138908447, -1.38780282263764,
              -1.387805349423,   -1.38882047800373, -1.49793700616569, -6.19975747800726,
              -22.197881249831,  -38.2046765807324, -38.2046769835162, -70.2261446634383,
              -54.2154105303641, -78.2315117307532, -86.2368787980699, -118.258347067337,
              -126.263714134653, -134.26908120197,  -158.28518240392,  -158.28518240392 } };


const double log_1_4 = std::log(1. / 4);
// based on Bryc, K., Patterson, N., and Reich, D. (2013). A novel approach
// to estimating heterozygosity from low-coverage genome sequence.
const double hetero_prior = 0.0005;
const double mut_prior = 1e-6;
const double homo_prior = 1 - hetero_prior - mut_prior; // not hetero and no somatic mutation
const double log_homo_prior = std::log(hetero_prior);


Filter::Filter(double theta, uint8_t cell_proportion)
    : theta(theta),
      cell_proportion(cell_proportion),
      log_theta_3(std::log(theta / 3)),
      log_one_minus_theta(std::log(1 - theta)) {
    log_factorial.reserve(171);
    log_factorial.push_back(1);
    for (uint32_t i = 1; i < 171; ++i) {
        log_factorial.push_back(log_factorial.back() * i);
    }
    for (uint32_t i = 0; i < 171; ++i) {
        log_factorial[i] = std::log(log_factorial[i]);
    }
}

double Filter::log_fact(uint32_t n) {
    return (n > 170) ? 0.5 * std::log(2 * M_PI * n) + n * std::log(n / M_E) : log_factorial.at(n);
}

double round_nearest_even(double x) {
    std::fesetround(FE_TONEAREST);
    return std::nearbyint(x);
}

bool Filter::is_two_sigmas_away(uint32_t coverage, std::array<uint16_t, 4> &base_count) {
    return base_count[3] > 2 * base_count[2] && base_count[2] > 2
            && base_count[2] > (coverage * theta) + 3 * sqrt(coverage * theta * (1 - theta));
}


bool Filter::is_significant(std::array<uint16_t, 4> &base_count) {
    // total coverage
    uint32_t coverage = base_count[0] + base_count[1] + base_count[2] + base_count[3];

    // if there are no reads or if we have just one base, do not keep
    if (coverage < 2) {
        return false;
    }

    // sort base_count; sorts in ascending order
    std::sort(base_count.begin(), base_count.end());

    if (base_count[2] == 0) { // all bases are the same
        return false;
    }

    // as simple as it looks, this condition improves results on the synthetic dataset
    // see
    // https://github.com/ratschlab/secedo-evaluation/commit/9822fcc63e028011eaa1e60c398101021e0eb64c
    if (base_count[2] + base_count[1] + base_count[0] < 5) {
        return false;
    }

    if (base_count[3] < 1.5 * base_count[2]) {
        return false;
    }

    // choose K threshold for the closest coverage, rounding to nearest even to emulate Python
    uint32_t threshold_idx
            = static_cast<uint32_t>(std::clamp(round_nearest_even((coverage / 10.)) - 1, 0., 19.));

    // log P(base_count|wt) for wt=homo
    double log_prob_homozygous
            = base_count[3] * log_one_minus_theta + (coverage - base_count[3]) * log_theta_3;
    // add prior on the most probable genotype (1/4, because all four homozygous genotypes are
    // equally likely)
    log_prob_homozygous += log_1_4;
    // add prior on null hypothesis (0.998)
    log_prob_homozygous += log_homo_prior;
    // the normalizing coefficient (probability of the evidence P(c1,c2,c3,c4)
    // double log_evidence = log_fact(coverage + 3) - log_6 - log_fact(coverage);
    // 1. All true allelles are c1
    double prob_all_c1 = homo_prior * std::pow(1 - theta, base_count[3])
            * std::pow(theta / 3, coverage - base_count[3]);

    // 2. The locus is heterozygous c1 c2
    double prob_hetero = hetero_prior * std::pow(0.5 - theta / 3, base_count[3] + base_count[2])
            * std::pow(theta / 3, base_count[0] + base_count[1]);
    // 3. The locus is homzygous + somatic mutation
    double prob_homo_som = homo_prior * mut_prior * std::pow(0.75 - 2 * theta / 3, base_count[3])
            * std::pow(0.25, base_count[2]) * std::pow(theta / 3, base_count[0] + base_count[1]);

    // 4. The locus is heterozygous + somatic mutation
    double prob_hetero_som = hetero_prior * mut_prior * std::pow(0.5 - theta, base_count[3])
            * std::pow(0.25, base_count[1] + base_count[2]) * std::pow(theta / 3, base_count[0]);
    // 5. Two somatic mutations
    double prob_two_somatic = hetero_prior * mut_prior * mut_prior * std::pow(1 - theta, coverage);
    double log_evidence
            = log(prob_all_c1 + prob_hetero + prob_homo_som + prob_hetero_som + prob_two_somatic);
    return log_prob_homozygous - log_evidence < Ks[cell_proportion].at(threshold_idx);
}

bool Filter::is_significant(const PosData &pos_data, uint16_t *coverage) {
    std::array<uint16_t, 4> base_count = { 0, 0, 0, 0 };
    for (uint32_t i = 0; i < pos_data.size(); ++i) {
        base_count[pos_data.base(i)]++;
    }
    *coverage = base_count[0] + base_count[1] + base_count[2] + base_count[3];
    return is_significant(base_count);
}

std::pair<std::vector<std::vector<PosData>>, double>
Filter::filter(const std::vector<std::vector<PosData>> &pos_data,
               const std::vector<uint32_t> &id_to_pos,
               const std::string &marker,
               uint32_t num_threads) {
    std::vector<std::vector<PosData>> result(pos_data.size());
    // one per thread to avoid lock contention
    std::vector<uint32_t> coverage_chr(pos_data.size());
    // can be atomic, very little lock contention
    std::atomic<uint32_t> total_positions = 0;
    std::ignore = num_threads;
#pragma omp parallel for num_threads(num_threads)
    for (uint32_t chr_idx = 0; chr_idx < pos_data.size(); ++chr_idx) {
        std::vector<PosData> filtered_data;
        for (uint32_t pos_idx = 0; pos_idx < pos_data[chr_idx].size(); ++pos_idx) {
            std::vector<uint32_t> read_ids;
            std::vector<uint16_t> cell_ids_and_bases;
            const PosData &pd = pos_data[chr_idx][pos_idx];
            std::array<uint16_t, 4> base_count = { 0, 0, 0, 0 };
            for (uint32_t cell_idx = 0; cell_idx < pd.size(); ++cell_idx) {
                if (id_to_pos[pd.group_id(cell_idx)] == NO_POS) {
                    continue;
                }
                read_ids.push_back(pd.read_ids[cell_idx]);
                cell_ids_and_bases.push_back(pd.group_ids_bases[cell_idx]);
                base_count[pd.base(cell_idx)]++;
            }

            if (is_significant(base_count)) {
                coverage_chr[chr_idx] += read_ids.size();
                filtered_data.emplace_back(pd.position, std::move(read_ids),
                                           std::move(cell_ids_and_bases));
            }
        }
        total_positions.fetch_add(filtered_data.size());
        result[chr_idx] = std::move(filtered_data);
    }

    uint32_t total_coverage = std::accumulate(coverage_chr.begin(), coverage_chr.end(), 0);
    double avg_coverage
            = total_positions == 0 ? 0 : static_cast<double>(total_coverage) / total_positions;
    logger()->trace("Avg coverage for cluster {}: {}. Total positions: {}", marker, avg_coverage,
                    total_positions);
    return { result, avg_coverage };
}
