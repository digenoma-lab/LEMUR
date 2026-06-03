# MergeBedMethyl

Stream-merge [modkit](https://github.com/nanoporetech/modkit) **bedMethyl** files from phased haplotypes into a single TSV matrix. Built for cohorts where each sample has `*_hp1.bedmethyl` and `*_hp2.bedmethyl` pairs.

## Output format

Each row is one genomic position (`chr`, `pos`). Each sample contributes six columns:

`{id}.hap1_counts`, `{id}.hap2_counts`, `{id}.hap1_cov`, `{id}.hap2_cov`, `{id}.hap1_percentage`, `{id}.hap2_percentage`

| chr | pos | CHI01.hap1_counts | … | CHI01.hap2_percentage | CHI02.hap1_counts | … |
|-----|-----|-------------------|---|-------------------------|-------------------|---|
| chr1 | 1015144 | 3 | … | 100 | 2 | … |

- **counts** / **cov**: from bedMethyl columns `N_modified` and valid coverage.
- **percentage**: `percent_modified` (0–100) from bedMethyl.
- **Missing** or low coverage: `.` for that haplotype field.

## Requirements

- C++17 compiler (g++ ≥ 7, clang ≥ 5)
- CMake ≥ 3.14

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
# or
make
```

Binary: `build/merge_bedmethyl`

```bash
cmake --install build   # optional, installs to CMAKE_INSTALL_PREFIX/bin
```

## Usage

```bash
merge_bedmethyl [-c N] [-s M] <output.tsv> <label1> <hp1> <hp2> [<label2> ...]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-c`, `--min-cov N` | Report methylation only if coverage (column 10) **>** N | `3` |
| `-s`, `--min-samples M` | Minimum samples with data per row | `N-1` (N = number of pairs) |
| `-h`, `--help` | Show help | |

### Example

```bash
./build/merge_bedmethyl merged.tsv \
  CHI01 results/modkit/CHI01.bed/CHI01_hp1.bedmethyl results/modkit/CHI01.bed/CHI01_hp2.bedmethyl \
  CHI02 results/modkit/CHI02.bed/CHI02_hp1.bedmethyl results/modkit/CHI02.bed/CHI02_hp2.bedmethyl
```

Require all samples at each site:

```bash
./build/merge_bedmethyl -s 2 -c 3 merged.tsv CHI01 ... CHI02 ...
```

## Slurm

```bash
export SAMPLES="CHI01 CHI02 CHI03"
export MODKIT_DIR=/path/to/results/modkit
export OUTPUT=merged.tsv
sbatch scripts/run_slurm.sh
```

## Tests

```bash
make test
# or
cd build && ctest --output-on-failure
```

## Project layout

```
MergeBedMethyl/
├── CMakeLists.txt
├── Makefile              # wrapper around CMake
├── include/merge_bedmethyl/
├── src/
├── tests/data/           # small bedMethyl fixtures
└── scripts/run_slurm.sh
```

## Input assumptions

- bedMethyl files are **sorted** by chromosome and start (modkit default).
- Column indices (0-based): start = 1, coverage = 9, percent modified = 10, N_modified = 11.

## License

MIT — see [LICENSE](LICENSE).
