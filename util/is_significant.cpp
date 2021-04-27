#include "is_significant.hpp"

#include <cfenv>

// thresholds K to use for pre-processing; values for coverage 10, 20, 30, ... 200
// always the largest values for any tumour proportion, rounded up to one decimal place
std::vector<double> Ks = { 0.872,   0.872,   0.872,   -1.238,  -1.358,   -5.749,  -10.139,
                           -13.912, -25.998, -33.936, -34.469, -54.084,  -44.453, -56.842,
                           -63.352, -88.297, -90.711, -96.841, -115.730, -112.601 };

std::vector<double> log_factorial;

void init() {
    log_factorial.reserve(171);
    log_factorial.push_back(1);
    for (uint32_t i = 1; i < 171; ++i) {
        log_factorial.push_back(log_factorial.back() * i);
    }
    for (uint32_t i = 0; i < 171; ++i) {
        log_factorial[i] = std::log(log_factorial[i]);
    }
}

double log_fact(uint32_t n) {
    static bool is_init = false;
    if (!is_init) {
        init();
        is_init = true;
    }
    return (n > 170) ? 0.5 * std::log(2 * M_PI * n) + n * std::log(n / M_E) : log_factorial.at(n);
}

double round_nearest_even(double x) {
    std::fesetround(FE_TONEAREST);
    return std::nearbyint(x);
}

bool is_significant(std::array<uint16_t, 4> &base_count, double theta) {
    double log_theta = std::log(theta / 3);
    double log_one_minus_theta = std::log(1 - theta);

    // sort base_count; sorts in ascending order
    std::sort(base_count.begin(), base_count.end());
    // total coverage
    uint32_t coverage = base_count[0] + base_count[1] + base_count[2] + base_count[3];

    // choose K for the closest coverage, rounding to nearest even to emulate Python
    uint32_t i = std::clamp(round_nearest_even((coverage / 10.) ) - 1, 0., 19.);

    // if there are no reads or if we have just one base, do not keep
    if (coverage == 0 || base_count[2] == 0) {
        return false;
    }
    // the multinomial coefficient
    double log_multinomial_coeff = log_fact(coverage) - log_fact(base_count[0])
            - log_fact(base_count[1]) - log_fact(base_count[2]) - log_fact(base_count[3]);
    // log P(base_count|wt) for wt=homo
    double log_prob_homozygous = log_multinomial_coeff + base_count[3] * log_one_minus_theta
            + (coverage - base_count[3]) * log_theta;
    // add prior on the most probable genotype (1/4, because all four homozygous genotypes are
    // equally likely)
    log_prob_homozygous += std::log(1. / 4);
    // add prior on null hypothesis (0.998)
    log_prob_homozygous += std::log(0.998);
    // the normalizing coefficient (normalizing for read depth)
    double log_normalizing_coef = log_fact(coverage + 3) - std::log(6) - log_fact(coverage);
    return log_normalizing_coef + log_prob_homozygous < Ks.at(i);
}

bool is_significant(const PosData &pos_data, double theta, uint32_t *coverage) {
    std::array<uint16_t, 4> base_count = { 0, 0, 0, 0 };
    for (const auto &cd : pos_data.cells_data) {
        base_count[cd.base]++;
    }
    *coverage = base_count[0] + base_count[1] + base_count[2] + base_count[3];
    return is_significant(base_count, theta);
}
