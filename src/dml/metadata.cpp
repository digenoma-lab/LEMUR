#include "dml/metadata.hpp"

#include "dml/util.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace dml {

namespace {

enum class DesignMode { Minimal, Ancestry, Pc };

std::string trim(const std::string& s) {
    std::size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) ++start;
    std::size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) --end;
    return s.substr(start, end - start);
}

int phenotype_bin_from_field(const std::string& phen, const std::string& case_label,
                             const std::string& control_label) {
    if (phen == case_label || phen == "2") return 1;
    if (phen == control_label || phen == "1") return 0;
    throw std::runtime_error("Unexpected phenotype value: " + phen);
}

DesignMode detect_design_mode(const std::vector<std::string>& header) {
    if (column_index(header, "PC1") >= 0) return DesignMode::Pc;
    if (column_index(header, "AMR") >= 0 && column_index(header, "COVERAGE_MEAN") >= 0) {
        return DesignMode::Ancestry;
    }
    return DesignMode::Minimal;
}

}  // namespace

CohortMeta load_metadata(const std::string& path, const std::vector<std::string>& meth_header,
                         bool sample_mode, const std::string& case_label,
                         const std::string& control_label) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open metadata: " + path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty metadata: " + path);

    const std::vector<std::string> header = split_csv_line(line, ',');
    const int id_idx = column_index(header, "sample_id");
    const int phen_idx = column_index(header, "phenotype");
    const int age_idx = column_index(header, "AGE");
    const int sex_idx = column_index(header, "SEX");
    const int bmi_idx = column_index(header, "BMI");
    const int cov_mean_idx = column_index(header, "COVERAGE_MEAN");
    const int amr_idx = column_index(header, "AMR");

    if (id_idx < 0 || phen_idx < 0 || age_idx < 0 || bmi_idx < 0) {
        throw std::runtime_error("Metadata must contain sample_id,phenotype,AGE,BMI columns");
    }

    const DesignMode mode = detect_design_mode(header);

    std::vector<int> pc_idx;
    if (mode == DesignMode::Pc) {
        if (sex_idx < 0 || cov_mean_idx < 0) {
            throw std::runtime_error("PC design requires SEX and COVERAGE_MEAN columns");
        }
        for (int k = 1; k <= 3; ++k) {
            const int idx = column_index(header, "PC" + std::to_string(k));
            if (idx < 0) {
                throw std::runtime_error("Metadata with PC1 must also contain PC2 and PC3");
            }
            pc_idx.push_back(idx);
        }
    } else if (mode == DesignMode::Ancestry) {
        if (sex_idx < 0 || cov_mean_idx < 0) {
            throw std::runtime_error("Ancestry design requires SEX and COVERAGE_MEAN columns");
        }
    }

    CohortMeta meta;
    meta.sample_mode = sample_mode;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto fields = split_csv_line(line, ',');
        if (static_cast<int>(fields.size()) <= bmi_idx) continue;

        SampleMeta s;
        s.sample_id = trim(fields[static_cast<std::size_t>(id_idx)]);
        const std::string phen = trim(fields[static_cast<std::size_t>(phen_idx)]);
        s.phenotype_bin = phenotype_bin_from_field(phen, case_label, control_label);

        const double age = std::stod(fields[static_cast<std::size_t>(age_idx)]);
        const double bmi = std::stod(fields[static_cast<std::size_t>(bmi_idx)]);

        if (mode == DesignMode::Minimal) {
            s.xrow = {1.0, static_cast<double>(s.phenotype_bin), age, bmi};
        } else if (mode == DesignMode::Ancestry) {
            const double sex = std::stod(fields[static_cast<std::size_t>(sex_idx)]);
            const double coverage = std::stod(fields[static_cast<std::size_t>(cov_mean_idx)]);
            s.xrow = {1.0, static_cast<double>(s.phenotype_bin), age, sex, bmi,
                      std::stod(fields[static_cast<std::size_t>(amr_idx)]), coverage};
        } else {
            const double sex = std::stod(fields[static_cast<std::size_t>(sex_idx)]);
            const double coverage = std::stod(fields[static_cast<std::size_t>(cov_mean_idx)]);
            s.xrow = {1.0, static_cast<double>(s.phenotype_bin), age, sex, bmi};
            for (int idx : pc_idx) {
                s.xrow.push_back(std::stod(fields[static_cast<std::size_t>(idx)]));
            }
            s.xrow.push_back(coverage);
        }

        std::string count_col;
        std::string cov_col;
        if (sample_mode) {
            count_col = s.sample_id + ".counts";
            cov_col = s.sample_id + ".cov";
        } else {
            throw std::runtime_error(
                "dml haplotype mode is not supported; use --sample with "
                "{id}.counts / {id}.cov columns");
        }
        s.count_idx = column_index(meth_header, count_col);
        s.cov_idx = column_index(meth_header, cov_col);
        if (s.count_idx < 0 || s.cov_idx < 0) {
            throw std::runtime_error("Methylation TSV missing columns for sample " + s.sample_id +
                                     " (expected " + count_col + ", " + cov_col + ")");
        }

        meta.samples.push_back(std::move(s));
    }

    if (meta.samples.empty()) {
        throw std::runtime_error("No samples loaded from metadata");
    }
    meta.n_params = static_cast<int>(meta.samples.front().xrow.size());
    return meta;
}

}  // namespace dml
