#include "variant_calling.hpp"

uint8_t likely_homozygous(const std::array<uint16_t, 4> &n_bases, double theta) {
    uint32_t cov = sum(n_bases.begin(), n_bases.end());
    if (cov < 9) {
        return NO_GENOTYPE;
    }
    auto max = std::max_element(n_bases.begin(), n_bases.end());
    // first term is expected bases, second term is standard deviation
    if (cov - *max <= std::round(cov * theta + std::sqrt(cov * theta * (1 - theta)))) {
        uint32_t base = max - n_bases.begin();
        return base + (base << 3);
    }
    return NO_GENOTYPE;
}

uint8_t most_likely_genotype(const std::array<uint16_t, 4> &n_bases,
                             const std::array<uint16_t, 4> &n_bases_total,
                             const std::array<uint32_t, 4> &nbases_total_idx,
                             double hetero_prior,
                             double theta) {
    // coverage for the given position; // we require coverage at least 9
    uint32_t cov = sum(n_bases.begin(), n_bases.end());

    double logTheta = std::log(theta / 3);
    double logOneMinusTheta = std::log(1 - theta);
    double logHalfMinusTheta = std::log(0.5 - theta / 3);

    std::vector<uint32_t> idx = argsort(n_bases.begin(), n_bases.end());

    if (cov < 9) {
        if (n_bases[idx[3]] >= cov - 1) {
            return (idx[3] << 3) | idx[3];
        }
        uint8_t b1 = idx[3];
        uint8_t b2 = idx[2];
        uint8_t b3 = idx[1];
        uint8_t bt1 = nbases_total_idx[3];
        uint8_t bt2 = nbases_total_idx[2];
        if (n_bases[b2] > n_bases[b3] && (b1 == bt1 || b1 == bt2) && (b2 == bt1 || b2 == bt2)) {
            return (idx[3] << 3) | idx[2];
        }
        return NO_GENOTYPE;
    }


    // probability of the most likely homozygous genotype
    double logProb_homo = n_bases[idx[3]] * logOneMinusTheta + (cov - n_bases[idx[3]]) * logTheta;
    // probability of the most likely heterozygous genotype
    double logProb_hetero = (n_bases[idx[2]] + n_bases[idx[3]]) * logHalfMinusTheta
            + (n_bases[idx[0]] + n_bases[idx[1]]) * logTheta + std::log(hetero_prior);

    // if logProb_homo == logProb_hetero, we are not able to decide
    if (logProb_homo == logProb_hetero) {
        return NO_GENOTYPE;
    }

    // homozygous genotype
    if (logProb_homo > logProb_hetero) {
        // check that the genotype is unique
        if (n_bases[idx[2]] == n_bases[idx[3]]) {
            return NO_GENOTYPE;
        }
        return (idx[3] << 3) | idx[3];
    }

    // heterozygous genotype
    if (n_bases[idx[2]] == n_bases[idx[1]]) // check uniqueness
        return NO_GENOTYPE;
    return (idx[3] << 3) | idx[2];
}

std::unordered_map<std::string, std::vector<ChrMap>> read_map(const std::string &map_file) {
    if (map_file.empty()) {
        return {};
    }
    if (!std::filesystem::exists(map_file)) {
        logger()->error("Map file: {} does not exist.", map_file);
        exit(1);
    }
    std::unordered_map<std::string, std::vector<ChrMap>> result;
    std::ifstream f(map_file);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::vector<std::string> cols = split(line, '\t');
        if (cols.size() != 8) {
            logger()->error("Invalid map file: {}. Has {} columns, expected 8.", map_file,
                            cols.size());
            std::exit(1);
        }
        if (cols[1][0] < '1' && cols[1][0] > '9' && cols[1][0] != 'X' && cols[1][0] != 'Y') {
            continue; // not a chromosome mapping
        }
        if (cols[6] == "SEQ") {
            continue;
        }
        uint32_t start_pos = std::stoul(cols[2]) - 1;
        uint32_t len = std::stoul(cols[0]);
        result[cols[1]].push_back(
                { chromosome_to_id(cols[3]), start_pos, len, cols[6] == "INS" ? 'I' : 'D' });
    }
    return result;
}

void apply_map(const std::vector<ChrMap> &map,
               const std::vector<uint8_t> &chr_data,
               std::vector<uint8_t> *new_chr_data) {
    uint32_t i = 0;
    new_chr_data->resize(0);
    for (const ChrMap &map_entry : map) {
        while (i < map_entry.start_pos) {
            new_chr_data->push_back(chr_data[i]);
            i++;
        }
        if (map_entry.tr == 'D') {
            for (uint32_t pos = 0; pos < map_entry.len; ++pos) {
                new_chr_data->push_back(CharToInt[(uint32_t)'N']);
            }
        } else {
            for (uint32_t pos = 0; pos < map_entry.len; ++pos) {
                i++;
            }
        }
    }
    while (i < chr_data.size()) {
        new_chr_data->push_back(chr_data[i]);
        i++;
    }
}

void read_contig(std::ifstream &fasta_file, std::vector<uint8_t> *chr_data) {
    std::string line;
    while (std::getline(fasta_file, line)) {
        for (char c : line) {
            chr_data->push_back(CharToInt[(int)c]);
        }
        if (fasta_file.peek() == '>') {
            break;
        }
    }
}

void get_next_chromosome(std::ifstream &fasta_file,
                         const std::unordered_map<std::string, std::vector<ChrMap>> &map,
                         bool is_diploid,
                         std::vector<uint8_t> *chr_data,
                         std::vector<uint8_t> *tmp1,
                         std::vector<uint8_t> *tmp2) {
    std::string line;
    if (!std::getline(fasta_file, line)) {
        return;
    }
    assert(line[0] == '>');
    char chromosome = line[1];
    std::string chr_name = line.substr(1);

    logger()->trace("Reading chromosome: {} (first allele)", chr_name);

    tmp1->resize(0);
    read_contig(fasta_file, tmp1);

    if (map.find(chr_name) != map.end()) {
        apply_map(map.at(chr_name), *tmp1, chr_data);
    } else {
        std::swap(*tmp1, *chr_data);
    }

    char ch1 = fasta_file.get();
    assert(ch1 == '>' || ch1 == EOF);

    char ch2 = fasta_file.peek();

    fasta_file.putback(ch1);

    // Varsim fasta chromosomes are called >X_maternal, >Y_paternal, etc.
    if (!is_diploid || (chromosome == 'X' && ch2 == 'Y') || chromosome == 'Y') {
        if (ch2 == 'Y') {
            logger()->info("Male genome detected; Simulating homozygous X");
        } else {
            logger()->info("Simulating homozygous {}", line);
        }
        for (uint32_t i = 0; i < chr_data->size(); ++i) {
            chr_data->at(i) |= (chr_data->at(i) << 3);
        }
        return;
    }

    std::getline(fasta_file, line); // skip the >chr_name line

    chr_name = line.substr(1);

    logger()->trace("Reading chromosome: {} (second allele)", line);
    tmp1->resize(0);
    read_contig(fasta_file, tmp1);

    if (map.find(chr_name) != map.end()) {
        apply_map(map.at(chr_name), *tmp1, tmp2);
    } else {
        std::swap(*tmp1, *tmp2);
    }

    if (chr_data->size() != tmp2->size()) {
        logger()->error(
                "Invalid reference genome. Maternal and paternal chromosome sizes don't match ({} "
                "vs {})",
                chr_data->size(), tmp2->size());
        std::exit(1);
    }

    for (uint32_t i = 0; i < chr_data->size(); ++i) {
        chr_data->at(i) = (chr_data->at(i) << 3) | tmp2->at(i);
    }
}

bool is_homozygous(uint8_t base) {
    return (base & 7) == (base >> 3);
}

void write_vcf_preamble(std::ofstream &f, const std::string &reference_file, uint32_t cluster) {
    f << "##fileformat=VCFv4.2" << std::endl;
    std::time_t end_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    f << "##fileDate=" << std::ctime(&end_time);
    f << "##source=SVC (Somatic Variant Caller)" << std::endl;
    f << "##reference=" << reference_file << std::endl;
    f << "##cluster=" << cluster << std::endl;
    f << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcluster-" + std::to_string(cluster)
      << std::endl;
}

/** Retturns true if a and b have the same genotype (irrespective of allele ordering) */
bool is_same_genotype(uint8_t a, uint8_t b) {
    return a == b || ((a >> 3) | ((a & 7) << 3)) == b;
}

/**
 * Get the bases that are not the same in the given genotypes. E.g. if the first is A/C and the
 * second is A/T,the result will be {C,T}.
 */
std::vector<std::pair<char, char>> get_differing_bases(uint8_t reference_genotype,
                                                       uint8_t genotype) {
    assert((reference_genotype & 7) < 6 && (reference_genotype >> 3) < 6 && (genotype & 7) < 6
           && (genotype >> 3) < 6);
    char r1 = IntToChar[reference_genotype & 7];
    char r2 = IntToChar[reference_genotype >> 3];
    if (r1 > r2) {
        std::swap(r1, r2);
    }

    char g1 = IntToChar[genotype & 7];
    char g2 = IntToChar[genotype >> 3];


    if (g1 > g2) {
        std::swap(g1, g2);
    }

    if (g1 == r1 && g2 == r2) {
        return {};
    }

    if (g1 == r1) {
        return { { r2, g2 } };
    } else if (g2 == r2) {
        return { { r1, g1 } };
    } else {
        return { { r1, g1 }, { r2, g2 } };
    }
}

bool check_is_diploid(std::ifstream &f) {
    std::string s;
    f >> s;
    f.seekg(0);
    return s.find("maternal") != s.npos;
}

void write_vcf_line(uint8_t chr_idx,
                    uint32_t position,
                    uint8_t reference_genotype,
                    uint8_t genotype,
                    const std::array<uint16_t, 4> nbases,
                    const std::string &info_format,
                    std::ofstream &vcf_file) {
    // write position to file if different
    if (!is_same_genotype(genotype, reference_genotype) && genotype != NO_GENOTYPE) {
        if (is_homozygous(reference_genotype)) {
            std::string alt(1, IntToChar[genotype & 7]);
            std::string gt = "1/1"; // VCF genotype indicator
            if (!is_homozygous(genotype)) {
                alt += IntToChar[genotype >> 3];
                gt = "0/1";
            }
            assert((reference_genotype & 7) < 6);
            // TODO: remove the nbases[] - debug only
            vcf_file << id_to_chromosome(chr_idx) << '\t' << position << "\t.\t"
                     << IntToChar[reference_genotype & 7] << '\t' << alt << info_format << gt
                     << "\t" << nbases[0] << " " << nbases[1] << " " << nbases[2] << " "
                     << nbases[3] << " " << std::endl;
        } else {
            std::vector<std::pair<char, char>> bases
                    = get_differing_bases(reference_genotype, genotype);
            for (const auto &p : bases) {
                vcf_file << id_to_chromosome(chr_idx) << '\t' << position << "\t.\t" << p.first
                         << '\t' << p.second << info_format << "1/1"
                         << "\t" << nbases[0] << " " << nbases[1] << " " << nbases[2] << " "
                         << nbases[3] << " " << std::endl;
            }
        }
    }
}

void variant_calling(const std::vector<std::vector<PosData>> &pos_data,
                     const std::vector<uint16_t> &clusters,
                     const std::string &reference_genome_file,
                     const std::string &map_file,
                     double hetero_prior,
                     double theta,
                     const std::filesystem::path &out_dir) {
    if (clusters.empty()) {
        return;
    }

    std::filesystem::create_directories(out_dir);
    if (!std::filesystem::exists(reference_genome_file)) {
        logger()->error("Reference genome {} does not exist", reference_genome_file);
        exit(1);
    }

    std::unordered_map<std::string, std::vector<ChrMap>> map = read_map(map_file);

    // clusters are numbered starting with 1; 0 means "no cluster"
    uint32_t num_clusters = *std::max_element(clusters.begin(), clusters.end()) + 1;
    std::ifstream fasta_file(reference_genome_file);
    bool is_diploid = check_is_diploid(fasta_file);
    std::vector<uint8_t> reference_chromosome, tmp1, tmp2;
    if (std::filesystem::file_size(reference_genome_file) > 1e6) { // not a test setting
        reference_chromosome.reserve(250'000'000);
        tmp1.reserve(250'000'000);
        tmp2.reserve(250'000'000);
    }

    std::vector<std::ofstream> vcfs(num_clusters);
    for (uint32_t i = 0; i < vcfs.size(); ++i) {
        vcfs[i].open(std::filesystem::path(out_dir) / ("cluster_" + std::to_string(i) + ".vcf"));
        write_vcf_preamble(vcfs[i], reference_genome_file, i);
    }
    std::ofstream vcf_global(std::filesystem::path(out_dir) / "common.vcf");

    std::vector<double> cell_scores(clusters.size());
    std::vector<uint32_t> cell_loci(clusters.size());


    // the qual,filter,info,format fields are identical for all rows.
    std::string info_format = "\t.\t.\tVARIANT_OVERALL_TYPE=SNP\tGT\t";

    std::ofstream f(out_dir / "variant");
    for (uint32_t chr_idx = 0; chr_idx < pos_data.size(); ++chr_idx) {
        const std::vector<PosData> &chromosome_data = pos_data[chr_idx];
        get_next_chromosome(fasta_file, map, is_diploid, &reference_chromosome, &tmp1, &tmp2);
        logger()->trace("Calling variants on chromosome {}", chr_idx + 1);
        for (const PosData &pd : chromosome_data) {
            // this happens when the chromosome has an insert right at the end
            if (pd.position - 1 >= reference_chromosome.size()) {
                break;
            }

            std::array<uint16_t, 4> n_bases_total = { 0, 0, 0, 0 };
            std::vector<std::array<uint16_t, 4>> nbases(num_clusters);
            for (uint32_t i = 0; i < pd.size(); ++i) {
                uint32_t cl_idx = std::round(clusters[pd.cell_id(i)]);
                if (std::abs(static_cast<int64_t>(cl_idx - clusters[pd.cell_id(i)])) <= 0.05) {
                    nbases[cl_idx][pd.base(i)]++;
                }
                n_bases_total[pd.base(i)]++;
                cell_loci[pd.cell_id(i)]++;
            }

            // check if this is likely a homozygous (germline) locus
            uint8_t pooled_homo_genotype = likely_homozygous(n_bases_total, theta);

            // -1 because PosData is 1-based to emulate samtools et al
            uint8_t reference_genotype = reference_chromosome[pd.position - 1];

            if (pooled_homo_genotype != NO_GENOTYPE && pooled_homo_genotype != reference_genotype) {
                uint8_t new_base = pooled_homo_genotype & 7;
                uint8_t ref_base = (reference_genotype & 7) == new_base
                        ? (reference_genotype >> 3) & 7
                        : reference_genotype & 7;
                vcf_global << id_to_chromosome(chr_idx) << '\t' << pd.position << "\t.\t"
                           << IntToChar[ref_base] << '\t' << IntToChar[new_base] << info_format
                           << "1/1"
                           << "\t" << n_bases_total[0] << " " << n_bases_total[1] << " "
                           << n_bases_total[2] << " " << n_bases_total[3] << " " << std::endl;
            }
            std::array<uint32_t, 4> n_bases_total_idx = argsort<uint16_t, 4>(n_bases_total);
            std::vector<uint8_t> genotypes(num_clusters, NO_GENOTYPE);
            std::vector<uint8_t> coverage(num_clusters);
            for (uint32_t cl_idx = 0; cl_idx < num_clusters; ++cl_idx) {
                genotypes[cl_idx] = most_likely_genotype(nbases[cl_idx], n_bases_total,
                                                         n_bases_total_idx, hetero_prior, theta);
            }
            bool all_same = true;
            // starting at 2, because cluster0 means "no cluster", so we don't care what genotype
            // those cells have
            for (uint32_t cl_idx = 2; cl_idx < num_clusters; ++cl_idx) {
                if (!is_same_genotype(genotypes[cl_idx - 1], genotypes[cl_idx])) {
                    all_same = false;
                    break;
                }
            }
            if (all_same) {
                // write position to the common VCF file if different from reference genotype
                write_vcf_line(chr_idx, pd.position, reference_genotype, genotypes[0],
                               n_bases_total, info_format, vcf_global);
                continue;
            }
            for (uint32_t cl_idx = 0; cl_idx < num_clusters; ++cl_idx) {
                // skip if it's same as the pooled genotype
                if (pooled_homo_genotype != NO_GENOTYPE
                    && genotypes[cl_idx] == pooled_homo_genotype) {
                    continue;
                }

                // write position to file if different
                write_vcf_line(chr_idx, pd.position, reference_genotype, genotypes[cl_idx],
                               nbases[cl_idx], info_format, vcfs[cl_idx]);
            }

            // now let's see which cells don't match the inferred genotype and penalize their score
            for (uint32_t i = 0; i < pd.size(); ++i) {
                uint32_t cluster = std::round(clusters[pd.cell_id(i)]);
                uint8_t base = pd.base(i);
                if (genotypes[cluster] != NO_GENOTYPE && base != (genotypes[cluster] & 7)
                    && base != (genotypes[cluster] >> 3)
                    && sum(nbases[cluster].begin(), nbases[cluster].end()) > 9) {
                    cell_scores[pd.cell_id(i)]++;
                }
            }
        }
    }
    for (uint32_t i = 0; i < cell_scores.size(); ++i) {
        cell_scores[i] /= cell_loci[i];
    }
    write_vec(std::filesystem::path(out_dir) / "scores", cell_scores);
}
