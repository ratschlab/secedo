#include "pileup.hpp"

#include "sequenced_data.hpp"
#include "util/logger.hpp"
#include "util/util.hpp"

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
#include <progress_bar/progress_bar.hpp>


#include <array>
#include <atomic>
#include <filesystem>
#include <string>
#include <vector>

// the maximum theoretical limit for an insert
constexpr uint32_t MAX_INSERT_SIZE = 1'000;

// number of BAM records processed at a time
constexpr uint32_t CHUNK_SIZE = 1'000'000;

// max number of files we are allowed to open at a time while creating pileups
// if MAX_OPEN_FILES < input_files.size(), we will merge files in batches and then merge the
// result
constexpr uint32_t MAX_OPEN_FILES = 100;

// human chromosome lengths (for progress reporting only); last two are X and Y
constexpr uint32_t chromosome_lengths[]
        = { 248'956'422, 242'193'529, 198'295'559, 190'214'555, 181'538'259, 170'805'979,
            159'345'973, 145'138'636, 138'394'717, 133'797'422, 135'086'622, 133'275'309,
            114'364'328, 107'043'718, 101'991'189, 90'338'345,  83'257'441,  80'373'285,
            58'617'616,  64'444'167,  46'709'983,  50'818'468,  156'040'895, 57'227'415 };

struct CellData {
    uint32_t read_id;
    uint16_t cell_id_and_base;

    inline uint8_t base() { return cell_id_and_base & 3; }
    inline uint16_t cell_id() { return cell_id_and_base >> 2; }
};

/**
 * Reads data between start_pos and end_pos from the specified BAM reader and places the result in
 * #data.
 * @return true if there is no more data available in this reader, false otherwise
 */
bool read_bam_file(const uint16_t cell_id,
                   uint32_t chromosome,
                   uint32_t max_coverage,
                   uint32_t min_base_quality,
                   uint32_t min_map_quality,
                   uint32_t min_alignment_score,
                   uint32_t start_pos,
                   uint32_t end_pos,
                   BamTools::BamReader &reader,
                   std::ofstream &out_map,
                   std::vector<std::vector<CellData>> *data,
                   std::vector<std::atomic<uint16_t>> *data_size,
                   std::unordered_map<std::string, uint64_t> *read_name_to_id,
                   std::atomic<uint64_t> *last_read_id) {
    // iterate through all alignments, only keeping ones with high map quality
    BamTools::BamAlignment al;

    while (reader.GetNextAlignmentCore(al)) {
        if (static_cast<uint32_t>(al.RefID) != chromosome) {
            return true; // all data for the given chromosome was processed
        }
        if (static_cast<uint32_t>(al.Position) >= end_pos) {
            return false; // current chunk was processed
        }

        assert(al.IsProperPair());
        assert(al.IsPaired());
        assert(!al.IsFailedQC());
        if (al.Position < static_cast<int32_t>(start_pos)) {
            assert(static_cast<uint32_t>(al.GetEndPosition()) >= start_pos);
            continue;
        }

        al.BuildCharData();

        auto read_id_iter = read_name_to_id->find(al.Name);
        uint32_t read_id;
        if (read_id_iter == read_name_to_id->end()) {
            read_id = (*last_read_id)++;
            (*read_name_to_id)[al.Name] = read_id;
            out_map << al.Name << "\t" << read_id << std::endl;
        } else {
            read_id = read_id_iter->second;
        }
        uint32_t offset = 0; // if the CIGAR string contains inserts, we need to adjust the offset
        uint32_t del_offset = 0;
        uint32_t cigar_idx = 0;
        // skip soft/hard clips
        while (al.CigarData[cigar_idx].Type == 'H' || al.CigarData[cigar_idx].Type == 'S') {
            cigar_idx++;
        }
        // the last position in the current cigar chunk
        uint32_t cigar_end = al.CigarData[cigar_idx].Length;
        for (uint32_t i = 0; i + offset < al.AlignedBases.size(); ++i) {
            // check if we stepped outside the current CIGAR chunk
            while (i >= cigar_end) {
                cigar_idx++;
                assert(cigar_idx < al.CigarData.size());

                // if the aligned string has inserted bases relative to the reference, we simply
                // skip those bases by increasing the offset and move on to the next CIGAR chunk,
                // if there is one
                if (al.CigarData.at(cigar_idx).Type == 'I') {
                    offset += al.CigarData[cigar_idx].Length;
                    // on rare occasions, the last chunk in a cigar string is an I; if that's the
                    // case we are done with the current alignment
                    if (i + offset >= al.AlignedBases.size()) {
                        assert(cigar_idx == al.CigarData.size() - 1);
                        break;
                    }
                    continue;
                } else if (al.CigarData[cigar_idx].Type == 'D') {
                    // deleted bases don't show up in AlignedBases, but they do show up in Qualities
                    del_offset += al.CigarData[cigar_idx].Length;
                }
                cigar_end += al.CigarData[cigar_idx].Length;
            }

            if (i + offset >= al.AlignedBases.size()) {
                assert(cigar_idx == al.CigarData.size() - 1);
                break;
            }

            uint8_t base = CharToInt[(uint8_t)al.AlignedBases[i + offset]];
            // make sure we have a '-' on a deleted position
            assert(al.CigarData[cigar_idx].Type != 'D' || al.AlignedBases[i + offset] == '-');

            if (base == 5
                || static_cast<uint32_t>(al.Qualities[i + offset - del_offset] - 33U)
                        < min_base_quality) {
                continue;
            }
            // if mapping quality is too low, skip the read
            if (al.MapQuality < min_map_quality) {
                continue;
            }
            // if mapping score is too low, skip the read
            uint32_t alignment_score = 0;
            al.GetTag("AS", alignment_score);
            if (alignment_score < min_alignment_score) {
                continue;
            }

            assert(al.Position + i < end_pos + MAX_INSERT_SIZE);
            uint32_t pos = al.Position + i - start_pos;
            std::vector<CellData> &cell_datas = data->at(pos);
            uint16_t current_coverage = data_size->at(pos).fetch_add(1);
            if (current_coverage >= max_coverage) {
                continue; // this position has suspiciously high coverage, treating as noise
            }
            uint16_t cell_id_base = (cell_id << 2) | (base & 3);
            cell_datas[current_coverage] = { read_id, cell_id_base };
        }
    }

    return true; // input was exhausted
}

/**
 * Reads and merges data from the given input files.
 * @return true if there is no more data available in this reader, false otherwise
 */
bool read_bam_chunk(const std::vector<std::filesystem::path> &input_files,
                    uint16_t batch_start,
                    uint32_t chromosome_id,
                    uint32_t max_coverage,
                    uint32_t min_base_quality,
                    uint32_t min_map_quality,
                    uint32_t min_alignment_score,
                    uint32_t num_threads,
                    uint32_t start_pos,
                    uint32_t end_pos,
                    std::ofstream &out_map,
                    std::vector<std::vector<CellData>> *data,
                    std::vector<std::atomic<uint16_t>> *data_size,
                    std::vector<std::unordered_map<std::string, uint64_t>> *read_name_to_id,
                    std::atomic<uint64_t> *last_read_id) {
    bool is_done = true;
// std::ignore = num_threads;
#pragma omp parallel for num_threads(num_threads)
    for (uint32_t i = 0; i < input_files.size(); ++i) {
        auto reader = std::make_unique<BamTools::BamReader>();
        if (!reader->Open(input_files[i])) {
            logger()->error("Could not open {}", input_files[i]);
            std::exit(1);
        }

        reader->OpenIndex(input_files[i].string() + ".bai");
        if (!reader->Jump(chromosome_id, start_pos)) {
            logger()->warn("Could not jump to chromosome {}, position {} in file {}", chromosome_id,
                           start_pos, input_files[i]);
            continue;
        }
        is_done &= read_bam_file(batch_start + i, chromosome_id, max_coverage, min_base_quality,
                                 min_map_quality, min_alignment_score, start_pos, end_pos, *reader,
                                 out_map, data, data_size, &(read_name_to_id->at(i)), last_read_id);
    }

    return is_done;
}

std::vector<std::vector<std::filesystem::path>>
get_batches(const std::vector<std::filesystem::path> &input_files, uint32_t batch_size) {
    std::vector<std::vector<std::filesystem::path>> result;
    auto beg = input_files.begin();
    uint32_t batch_count = 0;
    while (beg != input_files.end()) {
        auto end = beg + std::min(batch_size, static_cast<uint32_t>(input_files.end() - beg));
        std::vector<std::filesystem::path> filenames(beg, end);
        result.push_back(filenames);
        batch_count++;
        beg = end;
    }
    return result;
}

PosData create_pos_data(uint32_t pos, std::vector<CellData> cell_data, uint32_t size) {
    std::vector<uint16_t> cell_ids_bases(size);
    std::vector<uint32_t> read_ids(size);
    for (uint32_t i = 0; i < size; ++i) {
        cell_ids_bases[i] = cell_data[i].cell_id_and_base;
        read_ids[i] = cell_data[i].read_id;
    }
    return PosData(pos, std::move(read_ids), std::move(cell_ids_bases));
}

std::vector<PosData> pileup_bams(const std::vector<std::filesystem::path> &bam_files,
                                 const std::filesystem::path &out_pileup,
                                 bool write_text_file,
                                 uint32_t chromosome_id,
                                 uint32_t max_coverage,
                                 uint32_t min_base_quality,
                                 uint32_t min_map_quality,
                                 uint32_t min_alignment_score,
                                 uint32_t num_threads,
                                 uint16_t min_different) {
    size_t total_size = (chromosome_lengths[chromosome_id] / CHUNK_SIZE) + 1;

    logger()->trace("Allocating data structures...");
    std::atomic<uint64_t> last_read_id = 0; // used to map read names to shorter integer values
    std::vector<std::unordered_map<std::string, uint64_t>> read_name_to_id(bam_files.size());

    std::vector<std::vector<std::filesystem::path>> file_batches
            = get_batches(bam_files, MAX_OPEN_FILES);

    std::ofstream out_text(out_pileup.string() + ".txt");
    std::ofstream out_bin(out_pileup.string() + ".bin", std::fstream::binary);
    std::ofstream out_map(out_pileup.string() + ".map", std::fstream::binary);

    std::vector<PosData> result; // the final result containing the merged data

    // this is where each thread will write the parsed data for the current chunk
    std::vector<std::vector<CellData>> data(CHUNK_SIZE + MAX_INSERT_SIZE,
                                            std::vector<CellData>(max_coverage));
    // number of elements in each data vector
    std::vector<std::atomic<uint16_t>> data_size(CHUNK_SIZE + MAX_INSERT_SIZE);
    std::for_each(data_size.begin(), data_size.end(), [](auto &v) { v = 0; });

    uint32_t pos_count = 0;
    uint32_t start_pos = 0;
    // traverse in chunks of size CHUNK_SIZE until all data was read
    bool is_done = false;
    ProgressBar read_progress(total_size, "Reading progress", std::cout);
    while (!is_done) {
        is_done = true;
        for (uint32_t i = 0; i < file_batches.size(); ++i) {
            is_done &= read_bam_chunk(file_batches[i], i * MAX_OPEN_FILES, chromosome_id,
                                      max_coverage, min_base_quality, min_map_quality,
                                      min_alignment_score, num_threads, start_pos,
                                      start_pos + CHUNK_SIZE, out_map, &data, &data_size,
                                      &read_name_to_id, &last_read_id);
        }

        for (uint32_t pos = 0; pos < CHUNK_SIZE; ++pos) {
            uint16_t coverage = data_size[pos];
            if (coverage < 2 || coverage >= max_coverage) {
                continue; // positions with too low or too high coverage are ignored
            }
            std::array<uint16_t, 4> nbases = { 0, 0, 0, 0 };
            for (uint32_t i = 0; i < coverage; ++i) {
                nbases[data[pos][i].base()]++;
            }
            uint16_t max_bases = *std::max_element(nbases.begin(), nbases.end());

            if (coverage - max_bases
                < min_different) { // assuming homozgous germline, thus irrelevant
                continue;
            }
            pos_count++;
            uint32_t chromosome = chromosome_id + 1;
            // adding 1 to pos to emulate samtools, which is 1-based
            uint32_t position = start_pos + pos + 1;


            PosData pd = create_pos_data(position, data[pos], coverage);
            result.push_back(pd);

            if (write_text_file) {
                out_text << chromosome << '\t' << position << '\t' << coverage << '\t';
                std::sort(data[pos].begin(), data[pos].begin() + coverage,
                          [](CellData &a, CellData &b) { return a.cell_id() < b.cell_id(); });

                for (uint32_t i = 0; i < coverage; ++i) {
                    out_text << IntToChar[data[pos][i].base()];
                }
                out_text << '\t';
                for (uint32_t i = 0; i < coverage - 1U; ++i) {
                    out_text << data[pos][i].cell_id() << ",";
                }
                out_text << data[pos][coverage - 1].cell_id();

                out_text << '\t';
                for (uint32_t i = 0; i < coverage - 1U; ++i) {
                    out_text << data[pos][i].read_id << ",";
                }
                out_text << data[pos][coverage - 1].read_id;

                out_text << std::endl;
            }
            out_bin.write(reinterpret_cast<char *>(&position), sizeof(position));
            out_bin.write(reinterpret_cast<char *>(&coverage), sizeof(coverage));
            out_bin.write(reinterpret_cast<char *>(pd.read_ids.data()),
                          coverage * sizeof(pd.read_ids[0]));
            out_bin.write(reinterpret_cast<char *>(pd.group_ids_bases.data()),
                          coverage * sizeof(pd.group_ids_bases[0]));
        }

        for (uint32_t i = 0; i < MAX_INSERT_SIZE; ++i) {
            data_size[i].store(data_size[i + CHUNK_SIZE]);
            std::swap(data[i], data[i + CHUNK_SIZE]);
        }
        for (uint32_t i = MAX_INSERT_SIZE; i < data.size(); ++i) {
            data_size[i] = 0;
        }
        start_pos += CHUNK_SIZE;
        read_progress += 1;
    }
    logger()->info("Written {} positions to {}.txt/.bin", pos_count, out_pileup);
    return result;
}
