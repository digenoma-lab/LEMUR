# MergeBedMethyl

Stream-merge [modkit](https://github.com/nanoporetech/modkit) **bedMethyl** files from phased haplotypes into a single TSV matrix. Built for cohorts where each sample has `*_hp1.bedmethyl` and `*_hp2.bedmethyl` pairs.

## Output format

Tab-separated values (TSV). The first two columns are shared across all samples:

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `pos` | Start position (bedMethyl start, 0-based) |

For each sample label `{id}` given on the command line, six columns are appended **in this order**:

| Column | Source (bedMethyl) | Description |
|--------|-------------------|-------------|
| `{id}.hap1_counts` | `N_modified` (haplotype 1) | Methylated read count |
| `{id}.hap2_counts` | `N_modified` (haplotype 2) | Methylated read count |
| `{id}.hap1_cov` | valid coverage (haplotype 1) | Coverage at the locus |
| `{id}.hap2_cov` | valid coverage (haplotype 2) | Coverage at the locus |
| `{id}.hap1_percentage` | `percent_modified` (haplotype 1) | Methylation 0–100 |
| `{id}.hap2_percentage` | `percent_modified` (haplotype 2) | Methylation 0–100 |

With two samples `S1` and `S2`, the header is:

```
chr	pos	S1.hap1_counts	S1.hap2_counts	S1.hap1_cov	S1.hap2_cov	S1.hap1_percentage	S1.hap2_percentage	S2.hap1_counts	S2.hap2_counts	S2.hap1_cov	S2.hap2_cov	S2.hap1_percentage	S2.hap2_percentage
```

Example data row (from the test fixtures, `-c 3 -s 2`):

```
chr1	100	2	5	5	5	50	100	3	3	4	6	75	50
```

Meaning at `chr1:100`:

- **S1**: hap1 → 2 methylated reads, cov 5, 50%; hap2 → 5 reads, cov 5, 100%
- **S2**: hap1 → 3 reads, cov 4, 75%; hap2 → 3 reads, cov 6, 50%

### Missing values

If a haplotype has no row at the locus, or its coverage does not pass `-c` (coverage must be **>** N), all six fields for that haplotype are written as `.` (e.g. `.\t.\t.\t.\t.\t.` for both haplotypes of one sample).

A row is emitted only when at least `-s` samples have data at that locus (at least one haplotype per sample above the coverage threshold).

### Percentage formatting

Values come from bedMethyl `percent_modified` (0–100). `0` and `100` are written without decimals; other values drop trailing zeros (e.g. `75`, `50`, `12.5`).

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
merge_bedmethyl [-c N] [-s M] <output.tsv> <label1> <hp1> <hp2> [<label2> <hp3> <hp4> ...]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-c`, `--min-cov N` | Include haplotype fields only if valid coverage (column 9) **>** N | `3` |
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

The `merge_two_samples` test writes `build/tests/test_out.tsv` with the header and row shown in [Output format](#output-format).

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
- modkit bedMethyl column indices used (0-based): start = 1, valid coverage = 9, `percent_modified` = 10, `N_modified` = 11.

## License

MIT — see [LICENSE](LICENSE).
