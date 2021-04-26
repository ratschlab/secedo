#include "bam_reader.hpp"

#include "heap.hpp"
#include "preprocess.hpp"
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

constexpr uint32_t chromosome_lengths[]
        = { 248'956'422, 242'193'529, 198'295'559, 190'214'555, 181'538'259, 170'805'979,
            159'345'973, 145'138'636, 138'394'717, 133'797'422, 135'086'622, 133'275'309,
            114'364'328, 101'991'189, 90'338'345,  83'257'441,  80'373'285,  58'617'616,
            64'444'167,  46'709'983,  50'818'468,  156'040'895, 57'227'415 };



/**
 * Reads data between start_pos and end_pos from the specified BAM reader and places the result in
 * #data.
 * @return true if there is no more data available in this reader, false otherwise
 */
bool read_bam_file(const uint16_t cell_id,
                   uint32_t chromosome,
                   uint32_t max_coverage,
                   uint32_t min_base_quality,
                   uint32_t start_pos,
                   uint32_t end_pos,
                   BamTools::BamReader &reader,
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
        uint64_t read_id = (read_id_iter == read_name_to_id->end()) ? (*last_read_id)++
                                                                    : read_id_iter->second;
        uint32_t offset = 0; // if the CIGAR string contains inserts, we need to adjust the offset
        uint32_t cigar_idx = 0;
        uint32_t cigar_end = al.CigarData[0].Length;
        for (uint32_t i = 0; i < al.AlignedBases.size() - offset; ++i) {
            if (i >= cigar_end) {
                cigar_idx++;
                cigar_end = i + al.CigarData[cigar_idx].Length;
            }
            if (al.CigarData[cigar_idx].Type == 'I') {
                offset++;
            }

            uint8_t base = CharToInt[(uint8_t)al.AlignedBases[i + offset]];
            // make sure we have a '-' on a deleted position
            assert(al.CigarData[cigar_idx].Type != 'D' || base == 5);

            if (base == 5 || static_cast<uint32_t>(al.Qualities[i] - 33U) < min_base_quality) {
                continue;
            }

            assert(al.Position + i < end_pos + MAX_INSERT_SIZE);
            uint32_t pos = al.Position + i - start_pos - offset;
            std::vector<CellData> &cell_datas = data->at(pos);
            uint16_t current_coverage = data_size->at(pos).fetch_add(1);
            if (current_coverage >= max_coverage) {
                continue; // this position has suspiciously high coverage, treating as noise
            }
            // TODO: use integer ids
            cell_datas[current_coverage] = { std::to_string(read_id), cell_id, base };
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
                    uint32_t num_threads,
                    uint32_t start_pos,
                    uint32_t end_pos,
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
            logger()->error("Could not open input BAM files.");
            std::exit(1);
        }

        reader->OpenIndex(input_files[i].string() + ".bai");
        if (!reader->Jump(chromosome_id, start_pos)) {
            logger()->warn("Could not jump to chromosome {}, position {} in file {}", chromosome_id,
                           start_pos, input_files[i]);
            continue;
        }
        is_done &= read_bam_file(batch_start + i, chromosome_id, max_coverage, min_base_quality,
                                 start_pos, end_pos, *reader, data, data_size,
                                 &(read_name_to_id->at(i)), last_read_id);
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

std::vector<PosData> read_bam(const std::vector<std::filesystem::path> &input_files,
                              const std::filesystem::path &outfile,
                              uint32_t chromosome_id,
                              uint32_t max_coverage,
                              uint32_t min_base_quality,
                              uint32_t num_threads,
                              double sequencing_error_rate) {
    size_t total_size = (chromosome_lengths[chromosome_id] / CHUNK_SIZE) + 1;

    logger()->trace("Allocating data structures...");
    std::atomic<uint64_t> last_read_id = 0; // used to map read names to shorter integer values
    std::vector<std::unordered_map<std::string, uint64_t>> read_name_to_id(input_files.size());

    std::vector<std::vector<std::filesystem::path>> file_batches
            = get_batches(input_files, MAX_OPEN_FILES);

    std::ofstream fout(outfile);

    std::vector<PosData> result; // the final result containing the merged data

    // this is where each thread will write the parsed data for the current chunk
    std::vector<std::vector<CellData>> data(CHUNK_SIZE + MAX_INSERT_SIZE,
                                            std::vector<CellData>(max_coverage));
    // number of elements in each data vector
    std::vector<std::atomic<uint16_t>> data_size(CHUNK_SIZE + MAX_INSERT_SIZE);
    std::for_each(data_size.begin(), data_size.end(), [](auto &v) { v = 0; });

    uint64_t start_pos = 0;
    // traverse in chunks of size CHUNK_SIZE until all data was read
    bool is_done = false;
    ProgressBar read_progress(total_size, "Reading progress", std::cout);
    while (!is_done) {
        is_done = true;
        for (uint32_t i = 0; i < file_batches.size(); ++i) {
            is_done &= read_bam_chunk(file_batches[i], i * MAX_OPEN_FILES, chromosome_id,
                                      max_coverage, min_base_quality, num_threads, start_pos,
                                      start_pos + CHUNK_SIZE, &data, &data_size, &read_name_to_id,
                                      &last_read_id);
        }

        for (uint32_t pos = 0; pos < CHUNK_SIZE; ++pos) {
            if (data_size[pos] < 2 || data_size[pos] >= max_coverage) {
                continue; // positions with too low or too high coverage are ignored
            }
            // count the number of bases at position pos
            std::array<uint16_t, 4> base_count = { 0, 0, 0, 0 };
            for (uint32_t i = 0; i < data_size[pos]; ++i) {
                base_count[data[pos][i].base]++;
            }

            if (!is_significant(base_count, sequencing_error_rate)) {
                continue;
            }

            // adding 1 to start pos, bc. samtools is 1-based
            fout << (chromosome_id + 1) << '\t' << start_pos + pos + 1 << '\t' << data_size[pos]
                 << '\t';
            std::sort(data[pos].begin(), data[pos].begin() + data_size[pos],
                      [](auto &a, auto &b) { return a.cell_id < b.cell_id; });

            for (uint32_t i = 0; i < data_size[pos]; ++i) {
                fout << IntToChar[data[pos][i].base];
            }
            fout << '\t';
            for (uint32_t i = 0; i < data_size[pos] - 1U; ++i) {
                fout << data[pos][i].cell_id << ",";
            }
            fout << data[pos][data_size[pos] - 1].cell_id;
            fout << std::endl;

            result.push_back({ start_pos + pos, data[pos] });
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

    return result;
}
