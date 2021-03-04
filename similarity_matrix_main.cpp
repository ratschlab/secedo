/* Generate the distance matrix from the pre-processed mpileup file. */

#include "util.hpp"

#include <gflags/gflags.h>

#include <fstream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <iostream>
#include <math.h>
#include <numeric>

DEFINE_double(seq_error_rate, 0.001, "Sequencing errors rate, denoted by theta");
DEFINE_double(mutation_rate,
              0,
              "epsilon, estimated frequency of mutated loci in the pre-processed data set");
// estimate of how many positions are actually homozygous germline, were only included because of
// sequencing (or alignment!) errors
DEFINE_double(
        hzygous_prob,
        0,
        "The probability that a loci is homozygous, (not filtered correctly in the first step");

DEFINE_string(mpileup_file,
              "",
              "Input file containing 'pileup' textual format from an alignment, as written by "
              "preprocessing.py");
DEFINE_uint32(num_cells, 0, "Number of sequenced cells");
DEFINE_string(cells_file,
              "",
              "File with identifiers of cells (numbers between 0 and (num_of_cells - 1); the "
              "matrices will be computed only for these cells; if absent, all cells are used.");
DEFINE_string(cell_group_file,
              "",
              "For each cell, it contains a group to which it belongs. Cells from one group will "
              "be treated as one cell. Useful e.g. when creating data with artificially higher "
              "coverage; if absent, each cell has its own group");

DEFINE_uint32(run_id,
              0,
              "Identifier of the run (e.g. chromosome number); the resulting files will be called "
              "'mat_same_<run_id>.csv' and 'mat_diff_<id>.csv'");
DEFINE_uint32(max_insert, 1000, "Maximal considered insert size (for paired-end sequencing)");


/** Caches a bunch of combinations and powers used again and again in the computation */
struct Cache {
    const double epsilon = FLAGS_mutation_rate;
    const double h = FLAGS_hzygous_prob;

    // probability that two same letters will be read as different
    const double theta = FLAGS_seq_error_rate;
    const double theta2 = theta * theta;
    const double p_same_diff = 2 * theta * (1 - theta) + 2 * theta2 / 3;
    // probability that two same letters will be read as same
    const double p_same_same = 1 - p_same_diff;
    // probability that two different letters will be read as same
    const double p_diff_same = 2 * (1 - theta) * theta / 3 + 2 * theta2 / 9;
    // probability that two different letters will be read as different
    const double p_diff_diff = 1 - p_diff_same;

    std::vector<double> pow_p_same_same = { 1, p_same_same };
    std::vector<double> pow_p_same_diff = { 1, p_same_diff };
    std::vector<double> pow_p_diff_same = { 1, p_diff_same };
    std::vector<double> pow_p_diff_diff = { 1, p_diff_diff };

    std::vector<double> pow_1_epsilon_h = { 1 - epsilon - h };
    std::vector<double> pow_h = { 1, FLAGS_hzygous_prob };
    std::vector<double> pow_epsilon = { 1, epsilon };
    std::vector<double> pow_0_5 = { 1, 0.5 };
    std::vector<double> pow_pss_pds = { 1, p_same_same + p_diff_same };
    std::vector<double> pow_psd_pdd = { 1, p_same_diff + p_diff_diff };

    std::vector<std::vector<uint64_t>> comb = { { 1 }, { 1, 1 } };

    void update(uint32_t max) {
        if (max <= pow_p_same_same.size()) {
            return; // cache is up to date
        }
        for (uint32_t p = pow_p_same_same.size(); p < max; ++p) {
            pow_p_same_same.push_back(pow_p_same_same.back() * p_same_same);
            pow_p_same_diff.push_back(pow_p_same_diff.back() * p_same_diff);
            pow_p_diff_same.push_back(pow_p_diff_same.back() * p_diff_same);
            pow_p_diff_diff.push_back(pow_p_diff_diff.back() * p_diff_diff);
            pow_1_epsilon_h.push_back(pow_1_epsilon_h.back() * (1 - epsilon - h));
            pow_h.push_back(pow_h.back() * h);
            pow_epsilon.push_back(pow_epsilon.back() * epsilon);
            pow_0_5.push_back(pow_0_5.back() * 0.5);
            pow_pss_pds.push_back(pow_pss_pds.back() * (p_same_same + p_diff_same));
            pow_psd_pdd.push_back(pow_psd_pdd.back() * (p_same_diff + p_diff_diff));
            std::vector<uint64_t> new_comb(p + 2);
            new_comb[0] = 1;
            new_comb.back() = 1;
            for (uint32_t i = 1; i < comb.back().size(); ++i) {
                new_comb[i] = comb.back()[i - 1] + comb.back()[i];
            }
            comb.push_back(std::move(new_comb));
        }
    }
};

/**
 * Retrieve or compute the probability of result (x_s, x_d).
 */
double log_prob_of_result(uint32_t x_s,
                          uint32_t x_d,
                          double epsilon,
                          double h,
                          std::vector<std::vector<double>> &log_probs,
                          std::vector<std::vector<uint32_t>> &combs_xs_xd,
                          Cache &c) {
    if (log_probs[x_s][x_d] == std::numeric_limits<double>::max()) {
        c.update(x_s + x_d); // TODO: only update combinations to max
        double p = 0;
        for (uint32_t i = 0; i <= x_s; ++i) {
            for (uint32_t j = 0; j <= x_d; j++) {
                for (uint32_t k = 0; k <= x_s - i; ++k) {
                    for (uint32_t l = 0; l <= x_d - j; ++l) {
                        if (i + j > 0) {
                            p += c.comb[x_s][i] * c.comb[x_d][j] * c.comb[x_s - i][k]
                                    * c.comb[x_d - j][l] * c.pow_1_epsilon_h[i + j] * 0.5
                                    * (c.pow_p_same_same[i] * c.pow_p_same_diff[j]
                                       + c.pow_p_diff_same[i] * c.pow_p_diff_diff[j])
                                    * c.pow_h[k + l] * c.pow_p_same_same[k] * c.pow_p_same_diff[k]
                                    * c.pow_epsilon[x_s + x_d - i - j - k - l]
                                    * c.pow_0_5[x_s + x_d - i - j - k - l]
                                    * c.pow_pss_pds[x_s - i - k] * c.pow_psd_pdd[x_d - j - l];
                        } else {
                            p += c.comb[x_s][k] * c.comb[x_d][l] * c.pow_h[k + l]
                                    * c.pow_p_same_same[k] * c.pow_p_same_diff[k]
                                    * c.pow_epsilon[x_s + x_d - k - l]
                                    * c.pow_0_5[x_s + x_d - k - l] * c.pow_pss_pds[x_s - k]
                                    * c.pow_psd_pdd[x_d - l];
                        }
                    }
                }
            }
        }
        p = c.comb[x_s + x_d][x_s] * p;
        log_probs[x_s][x_d] = log(p);
    }
    combs_xs_xd[x_s][x_d] += 1;
    return log_probs[x_s][x_d];
}

void compare_with_reads(
        const std::unordered_map<uint32_t, std::tuple<std::string, uint32_t, uint32_t, uint32_t>>
                &active_reads,
        uint32_t r_id,
        const std::vector<uint32_t> &cell_ids,
        const std::vector<uint32_t> &cell_groups,
        std::vector<std::vector<double>> &mat_same,
        std::vector<std::vector<double>> &mat_diff,
        std::vector<std::vector<double>> &log_probs_same,
        std::vector<std::vector<double>> &log_probs_diff,
        std::vector<std::vector<uint32_t>> &combs_xs_xd,
        Cache &cache) {
    auto [seq, pos, cell_id, unused] = active_reads.at(r_id);
    const double epsilon = FLAGS_mutation_rate;
    const double h = FLAGS_hzygous_prob;

    for (const auto &[r_id_2, r_value_2] : active_reads) {
        if (r_id_2 == r_id) {
            std::cerr << "The same read id was present twice in active_reads; "
                         "something is wrong.";
            std::exit(1);
        }
        //  read is not from the same cell -> count the number of matches and mismatches
        //  in the overlap
        uint32_t cell_id_2 = std::get<2>(r_value_2);
        if (cell_id_2 != cell_id) {
            uint32_t pos2 = std::get<1>(r_value_2);
            std::string seq2 = std::get<0>(r_value_2);
            uint32_t x_s = 0;
            uint32_t x_d = 0;
            uint32_t diff = pos2 - pos;
            if (pos > pos2) {
                std::swap(pos, pos2);
                std::swap(seq, seq2);
                diff = pos - pos2;
            }

            for (uint32_t i = 0; i < seq2.size(); ++i) {
                if (seq2[i] != '*' && seq[i + diff] != '*'
                    && toupper(seq2[i]) == toupper(seq[i + diff])) {
                    x_s += 1; // both reads are the same
                } else if (seq2[i] != '*' && seq[i + diff] != '*') {
                    x_d += 1; // both reads are different
                }
            }

            // update the distance matrices
            if (x_s + x_d > 0) {
                // TODO: use a map instead of std::find
                uint32_t index1 = cell_groups[std::distance(
                        cell_ids.begin(), std::find(cell_ids.begin(), cell_ids.end(), cell_id))];
                uint32_t index2 = cell_groups[std::distance(
                        cell_ids.begin(), std::find(cell_ids.begin(), cell_ids.end(), cell_id_2))];
                if (index1 != index2) {
                    mat_same[index1][index2]
                            += log_prob_of_result(x_s, x_d, 0, h + 0.5 * epsilon, log_probs_same,
                                                  combs_xs_xd, cache);
                    mat_same[index2][index1] = mat_same[index1][index2];
                    mat_diff[index1][index2]
                            += log_prob_of_result(x_s, x_d, epsilon, h, log_probs_diff, combs_xs_xd,
                                                  cache);
                    mat_diff[index2][index1] = mat_diff[index1][index2];
                }
            }
        }
    }
}

/**
 * Compute mat_same (the matrix giving probabilities of cells i and j given they are in
 * the same cluster) and mat_diff (prob. of cells i and j given they are in different clusters)
 */
void computeSimilarityMatrix(const std::vector<uint32_t> &cell_ids,
                             const std::vector<uint32_t> &cell_groups) {
    constexpr uint32_t read_len = 100; // maximum read length for computing the
    Cache cache;

    // arrays with values of already computed probabilities, max_double if not yet computed
    std::vector<std::vector<double>> log_probs_same(
            2 * read_len, std::vector<double>(2 * read_len, std::numeric_limits<double>::max()));
    std::vector<std::vector<double>> log_probs_diff(
            2 * read_len, std::vector<double>(2 * read_len, std::numeric_limits<double>::max()));
    // count how many times we have seen a given combination of x_s, x_d
    std::vector<std::vector<uint32_t>> combs_xs_xd(2 * read_len,
                                                   std::vector<uint32_t>(2 * read_len, 0));

    //  a dictionary holding the active reads
    // key: read id
    // value: list with four values:
    //	- read sequence (string)
    //	- starting position (line number) (int)
    //	- cell id (string)
    //	- starting position in real coordinates (int)
    std::unordered_map<uint32_t, std::tuple<std::string, uint32_t, uint32_t, uint32_t>>
            active_reads;
    // distance matrices
    uint32_t num_cells = FLAGS_num_cells;
    std::vector<std::vector<double>> mat_same(num_cells, std::vector<double>(num_cells, 0));
    std::vector<std::vector<double>> mat_diff(num_cells, std::vector<double>(num_cells, 0));

    // line counter
    uint32_t n_lines = 0;
    // process the mpileup file line by line
    std::ifstream f(FLAGS_mpileup_file);
    std::string line;
    std::vector<uint32_t> line_cells;
    std::string line_bases;
    std::vector<uint32_t> line_reads;
    while (std::getline(f, line)) {
        n_lines++;
        std::stringstream line_stream(line);
        std::string s;
        uint32_t idx = 0;
        uint32_t line_pos = -1;
        // split the line by tabs
        while (getline(line_stream, s, '\t')) {
            switch (idx) {
                case 1: // 2nd field: position
                    line_pos = std::stoi(s);
                    break;
                case 3: // 4th field: bases
                    line_bases = s;
                    break;
                case 4: // 5th field: cell identifiers
                    line_cells = int_split(s, ',');
                    break;
                case 5: // 6th field: read identifiers
                    line_reads = int_split(s, ',');
                    break;
            };
            ++idx;
        }
        for (auto it = active_reads.begin(); it != active_reads.end();) {
            auto &[r_id, r_value] = *it;
            auto line_it = std::find(line_reads.begin(), line_reads.end(), r_id);
            if (line_it != line_reads.end()) { // the read continues
                uint32_t tmp_pos
                        = line_it - line_reads.begin(); // position of the read in line_reads
                char tmp_base = line_bases[tmp_pos]; // the sequenced base
                std::get<0>(r_value) += tmp_base;
                // if we are not more than max_insert from start of the mate,
                // it is possible the read still continues, only this particular base is unknown
                // (deletion, or low quality, or we are in between the two mates)
                ++it;
            } else if (std::get<3>(r_value) >= line_pos - FLAGS_max_insert) {
                std::get<0>(r_value) += '*';
                ++it;
            } else { // the read does not continue; compute its overlaps with all other active reads
                uint32_t cell_id = std::get<2>(r_value);
                uint32_t pos = std::get<1>(r_value);
                std::string seq = std::get<0>(r_value);
                // remove the read from active_reads
                it = active_reads.erase(it);

                // compare with all other reads in active_reads
                compare_with_reads(active_reads, r_id, cell_ids, cell_groups, mat_same, mat_diff,
                                   log_probs_same, log_probs_diff, combs_xs_xd, cache);
            }
        }

        // for each read in line_reads, check, if it is in active_reads; if not, add it
        for (uint32_t i = 0; i < line_reads.size(); ++i) {
            if (!active_reads.contains(line_reads[i])) {
                std::string b(1, line_bases[i]);
                uint32_t cell_id = line_cells[i];
                active_reads[line_reads[i]] = std::make_tuple(b, n_lines, cell_id, line_pos);
            }
        }
    }
    f.close();

    // in the end, process all the remaining active reads
    for (auto it = active_reads.begin(); it != active_reads.end();) {
        const auto &[r_id, r_value] = *it;
        it = active_reads.erase(it); // delete the read from active_reads

        // compare with all other reads in active_reads
        compare_with_reads(active_reads, r_id, cell_ids, cell_groups, mat_same, mat_diff,
                           log_probs_same, log_probs_diff, combs_xs_xd, cache);
    }

    // save mat_same and mat_diff into file
    write_mat("mat_same_" + std::to_string(FLAGS_run_id) + ".csv", mat_same);
    write_mat("mat_diff_" + std::to_string(FLAGS_run_id) + ".csv", mat_diff);
    write_mat("combs_xs_xd_" + std::to_string(FLAGS_run_id) + ".csv", combs_xs_xd);
}

int main(int argc, char *argv[]) {
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    std::vector<uint32_t> cells_ids;
    if (FLAGS_cells_file != "") {
        cells_ids = int_split(read_file(FLAGS_cells_file), ' ');
        // num_cells = len(cells_ids);
    } else {
        cells_ids.resize(FLAGS_num_cells);
        std::iota(cells_ids.begin(), cells_ids.end(), 0);
    }

    std::vector<uint32_t> cell_groups;
    if (FLAGS_cell_group_file != "") {
        cell_groups = int_split(read_file(FLAGS_cells_file), ',');
    } else {
        cells_ids.resize(FLAGS_num_cells);
        std::iota(cells_ids.begin(), cells_ids.end(), 0);
    }

    computeSimilarityMatrix(cells_ids, cell_groups);
}
