#include "bam_reader.hpp"

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

// the maximum theoretical limit for an insert
constexpr uint32_t MAX_INSERT_SIZE = 1'000;

// number of BAM records processed at a time
constexpr uint32_t CHUNK_SIZE = 100'000;

// max number of files we are allowed to open at a time while creating pileups
// if MAX_OPEN_FILES < input_files.size(), we will merge files in batches and then merge the
// result
constexpr uint32_t MAX_OPEN_FILES = 100;


/**
 * Reads data between start_pos and end_pos from the specified BAM reader and places the result in #data.
 * @return true if there is no more data available in this reader, false otherwise
 */
bool read_thread(const std::unordered_map<std::string, uint16_t> &fname_to_id,
                 uint32_t max_coverage,
                 uint32_t chromosome,
                 uint32_t start_pos,
                 uint32_t end_pos,
                 BamTools::BamMultiReader &reader,
                 std::vector<std::vector<CellData>> *data,
                 std::vector<std::atomic<uint16_t>> *data_size,
                 std::atomic<uint64_t> *last_read_id) {
    std::unordered_map<std::string, uint64_t> read_name_to_id;

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
        assert(al.Position >= static_cast<int32_t>(start_pos));

        // TODO: figure this out
        //        if (al.MapQuality < 90) {
        //            continue;
        //        }

        al.BuildCharData();
        auto read_id_iter = read_name_to_id.find(al.Name);
        uint64_t read_id;
        read_id = (read_id_iter == read_name_to_id.end()) ? (*last_read_id)++
                                                          : read_id_iter->second;
        for (uint32_t i = 0; i < al.AlignedBases.size(); ++i) {
            uint8_t base = CharToInt[(uint8_t)al.AlignedBases[i]];
            if (base == 5) {
                continue; // not ACGT, probably an inserted or deleted char
            }
            assert(al.Position + i < end_pos + MAX_INSERT_SIZE);
            std::vector<CellData> &cell_datas = data->at(al.Position + i - start_pos);
            uint16_t current_coverage = data_size->at(al.Position + i - start_pos).fetch_add(1);
            if (current_coverage >= max_coverage) {
                continue; // this position has suspiciously high coverage, treating as noise
            }
            // TODO: use integer ids
            cell_datas[current_coverage].cell_id = fname_to_id.at(al.Filename);
            cell_datas[current_coverage].read_id = std::to_string(read_id);
            cell_datas[current_coverage].base = base;
        }
    }

    return true; // input was exhausted
}

/**
 * Reads and merges data from the given input files.
 * @return the read and merged data
 */
std::vector<PosData> read_bam_chunk(const std::vector<std::filesystem::path> &input_files,
                                    const std::unordered_map<std::string, uint16_t> &fname_to_id,
                                    const std::filesystem::path &outfile,
                                    uint32_t chromosome_id,
                                    uint32_t max_coverage,
                                    uint32_t num_threads,
                                    double sequencing_error_rate) {
    uint32_t start_pos = 0;

    std::ofstream fout(outfile);

    std::vector<std::vector<CellData>> data(CHUNK_SIZE + MAX_INSERT_SIZE,
                                            std::vector<CellData>(max_coverage));
    // number of elements in each data vector
    std::vector<std::atomic<uint16_t>> data_size(CHUNK_SIZE + MAX_INSERT_SIZE);
    std::for_each(data_size.begin(), data_size.end(), [](auto &v) { v = 0; });

    std::atomic<uint64_t> last_read_id = 0; // used to map read names to shorter integer values

    size_t batch_size = std::max(1UL, input_files.size() / num_threads);

    logger()->trace("Reading chunks of ~{} files on {} threads", batch_size, num_threads);

    std::vector<std::unique_ptr<BamTools::BamMultiReader>> readers;

    size_t curr_pos = 0; // current position in the input files (for measuring progress)
    size_t end_pos = 0; // end position in the input files

    auto beg = input_files.begin();
    while (beg != input_files.end()) {
        auto end = beg + std::min(batch_size, static_cast<size_t>(input_files.end() - beg));
        std::vector<std::string> batch_fnames(beg, end);
        beg = end;

        auto reader = std::make_unique<BamTools::BamMultiReader>();
        if (!reader->Open(batch_fnames)) {
            logger()->error("Could not open input BAM files.");
            std::exit(1);
        }
        assert(reader->GetMergeOrder() == BamTools::BamMultiReader::MergeOrder::MergeByCoordinate);

        std::vector<std::string> index_fnames;
        for (const auto &fname : batch_fnames) {
            index_fnames.push_back(fname + ".bai");
        }
        reader->OpenIndexes(index_fnames);
        // jump to the beginning of next chromosome to figure out how much data we need to read
        reader->Jump(chromosome_id + 1);
        end_pos += reader->GetPosInFile();
        logger()->debug("End of next chromosome: {}", end_pos);
        // jump back to the current chromosome
        reader->Jump(chromosome_id);
        curr_pos += reader->GetPosInFile();
        logger()->debug("End of current chromosome: {}", curr_pos);
        readers.push_back(std::move(reader));
    }

    ProgressBar read_progress(end_pos - curr_pos, "Reading progress", std::cout);

    std::vector<PosData> result;
    while (true) {
        bool is_done = true;
#pragma omp parallel for num_threads(num_threads)
        for (uint32_t i = 0; i < readers.size(); ++i) {
            is_done &= read_thread(fname_to_id, max_coverage, chromosome_id, start_pos,
                                   start_pos + CHUNK_SIZE, *readers[i], &data, &data_size,
                                   &last_read_id);
        }

        for (uint32_t pos = 0; pos < CHUNK_SIZE; ++pos) {
            if (data_size[pos] < 2 || data_size[pos] >= max_coverage) {
                continue; // positions with too low or too high coverage are ignored
            }
            std::array<uint16_t, 4> base_count = { 0, 0, 0, 0 };
            for (uint32_t i = 0; i < data_size[pos]; ++i) {
                base_count[data[pos][i].base]++;
            }

            if (!is_significant(base_count, sequencing_error_rate)) {
                continue;
            }

            fout << "Size: " << data_size[pos] << "\t Start pos: " << start_pos + pos
                 << "\tCell ids: ";
            for (uint32_t i = 0; i < data_size[pos]; ++i) {
                fout << data[pos][i].cell_id << ",";
            }
            fout << "\tBases: ";

            for (uint32_t i = 0; i < data_size[pos]; ++i) {
                fout << IntToChar[data[pos][i].base];
            }
            fout << std::endl;

            result.push_back({ start_pos + pos, data[pos] });
        }

        if (is_done) {
            break;
        }

        for (uint32_t i = 0; i < MAX_INSERT_SIZE; ++i) {
            data_size[i].store(data_size[i + CHUNK_SIZE]);
            std::swap(data[i], data[i + CHUNK_SIZE]);
        }
        for (uint32_t i = MAX_INSERT_SIZE; i < data.size(); ++i) {
            data_size[i] = 0;
        }
        start_pos += CHUNK_SIZE;
        size_t new_file_pos = 0;
        for (const auto &r : readers) {
            new_file_pos += r->GetPosInFile();
        }
        read_progress += (new_file_pos - curr_pos);
        curr_pos = new_file_pos;
    }
    return result;
}

struct Cmp {
    int operator()(const PosData &a, const PosData &b) { return a.position < b.position; }
};

std::vector<PosData> merge(const std::vector<std::vector<PosData>> &sources) {
    MergeHeap<PosData, Cmp> heap;
    size_t total_size = 0;
    std::vector<size_t> source_idxs(sources.size(), 0);
    for (uint32_t i = 0; i < sources.size(); ++i) {
        total_size += sources[i].size();
        if (!sources[i].empty()) {
            heap.emplace(sources[i][0], i);
            source_idxs[i] = 1;
        }
    }
    std::vector<PosData> result;
    result.resize(total_size / 20); // because we expect a 0.05 coverage

    while (!heap.empty()) {
        auto [el, source_index] = heap.pop();

        if (result.empty() || result.back().position != el.position) {
            result.push_back(el);
        } else {
            result.back().cells_data.insert(result.back().cells_data.end(), el.cells_data.begin(),
                                            el.cells_data.end());
        }
        size_t cur_idx = ++source_idxs[source_index];
        if (cur_idx < sources.size()) {
            heap.emplace(sources[source_index][cur_idx], source_index);
        }
    }
    return result;
}

std::vector<PosData> read_bam(const std::vector<std::filesystem::path> &input_files,
                              const std::filesystem::path &outfile,
                              uint32_t chromosome_id,
                              uint32_t max_coverage,
                              uint32_t num_threads,
                              double sequencing_error_rate) {
    uint32_t cell_id = 0;
    std::unordered_map<std::string, uint16_t> fname_to_id;
    for (const auto &fname : input_files) {
        fname_to_id[fname] = cell_id++;
    }

    std::vector<std::vector<std::string>> file_batches;

    std::vector<std::vector<PosData>> pos_batches;

    auto beg = input_files.begin();
    uint32_t batch_count = 0;
    while (beg != input_files.end()) {
        logger()->info("Reading batch {}/{}", batch_count + 1,
                       std::max(1UL, input_files.size() / MAX_OPEN_FILES));
        auto end = beg + std::min(MAX_OPEN_FILES, static_cast<uint32_t>(input_files.end() - beg));
        std::vector<std::filesystem::path> filenames(beg, end);
        std::vector<PosData> pos_data
                = read_bam_chunk(filenames, fname_to_id, outfile, chromosome_id, max_coverage,
                                 num_threads, sequencing_error_rate);
        pos_batches.push_back(pos_data);
        batch_count++;
        beg = end;
    }

    logger()->info("Merging {} batches...", batch_count);
    return merge(pos_batches);
}
