#pragma once

#include "dml/metadata.hpp"
#include "dml/model.hpp"

#include <fstream>
#include <string>
#include <vector>

namespace dml {

class MethylationStream {
   public:
    MethylationStream(const std::string& path, const CohortMeta& meta);

    bool read_batch(std::vector<SiteInput>& batch, std::size_t max_rows);

   private:
    std::ifstream in_;
    const CohortMeta& meta_;
    std::vector<std::string> header_;
    bool header_read_ = false;
};

std::vector<std::string> read_tsv_header(const std::string& path);

}  // namespace dml
