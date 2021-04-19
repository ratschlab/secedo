#pragma once

#include "logger.hpp"
#include "preprocess.hpp"
#include "sequenced_data.hpp"
#include "util.hpp"

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>

#include <array>
#include <atomic>
#include <filesystem>
#include <string>
#include <vector>

uint32_t MAX_INSERT_SIZE = 1'000;

bool read_thread(const std::unordered_map<std::string, uint16_t> &fname_to_id,
                 uint32_t max_coverage,
                 uint32_t chromosome,
                 uint32_t start_pos,
                 uint32_t end_pos,
                 BamTools::BamMultiReader &reader,
                 std::vector<std::vector<CellData>> *data,
                 std::vector<std::atomic<uint16_t>> *data_size,
                 std::atomic<uint64_t> *last_read_id,
                 const std::function<void()> &progress) {
    std::unordered_map<std::string, uint64_t> read_name_to_id;

    // iterate through all alignments, only keeping ones with high map quality
    BamTools::BamAlignment al;
    while (reader.GetNextAlignmentCore(al)) {
        if (static_cast<uint32_t>(al.RefID) != chromosome) {
            return true; // all data for the given chromosome was processed
        }
        if (static_cast<uint32_t>(al.Position) >= end_pos) {
            // TODO: this may not be necessary
            if (reader.GetNextAlignmentCore(al)) { // check if we exhausted the content
                reader.Jump(al.RefID, al.Position); // unread the last read
                return false; // this chunk is finished
            }
            return true; // all done
        }
        progress();

        assert(al.IsProperPair());
        assert(al.IsPaired());
        assert(!al.IsFailedQC());
        assert(al.Position >= static_cast<int32_t>(start_pos));

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

std::vector<PosData> read_bams(const std::vector<std::filesystem::path> &inputFilenames,
                               const std::filesystem::path &outfile,
                               uint32_t chromosome_id,
                               uint32_t max_coverage,
                               uint32_t num_threads,
                               double sequencing_error_rate,
                               const std::function<void()> &progress) {
    uint32_t start_pos = 0;
    uint32_t CHUNK_SIZE = 100'000;

    std::ofstream fout(outfile);

    std::vector<std::vector<CellData>> data(CHUNK_SIZE + MAX_INSERT_SIZE,
                                            std::vector<CellData>(max_coverage));
    // number of elements in each data vector
    std::vector<std::atomic<uint16_t>> data_size(CHUNK_SIZE + MAX_INSERT_SIZE);
    for (auto &v : data_size) {
        v = 0;
    }
    std::atomic<uint64_t> last_read_id = 0; // used to map read names to shorter integer values

    size_t batch_size = std::max(1UL, inputFilenames.size() / num_threads);

    logger()->trace("Reading {} files on {} threads", batch_size, num_threads);

    uint32_t cell_id = 0;
    std::unordered_map<std::string, uint16_t> fname_to_id;
    for (const auto &fname : inputFilenames) {
        fname_to_id[fname] = cell_id++;
    }

    std::vector<std::unique_ptr<BamTools::BamMultiReader>> readers;

    auto beg = inputFilenames.begin();
    while (beg != inputFilenames.end()) {
        auto end = beg + std::min(batch_size, static_cast<size_t>(inputFilenames.end() - beg));
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
        reader->Jump(chromosome_id);
        readers.push_back(std::move(reader));
    }
    std::vector<PosData> result;
    while (true) {
        logger()->trace("Next chunk");
        bool is_done = true;
        // #pragma omp parallel for num_threads(num_threads)
        for (uint32_t i = 0; i < readers.size(); ++i) {
            is_done &= read_thread(fname_to_id, max_coverage, chromosome_id, start_pos,
                                   start_pos + CHUNK_SIZE, *readers[i], &data, &data_size,
                                   &last_read_id, progress);
        }

        for (uint32_t pos = 0; pos < CHUNK_SIZE; ++pos) {
            if (data_size[pos] < 2 || data_size[pos] >= max_coverage) {
                continue; // positions with too low or too high coverage are ignored
            }
            std::array<uint16_t, 4> base_count = { 0, 0, 0, 0 };
            for (uint32_t i = 0; i < data_size[pos]; ++i) {
                base_count[data[pos][i].base]++;
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

            if (!is_significant(base_count, sequencing_error_rate)) {
                continue;
            }

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
    }
    return result;
}
