
# LEMUR
### Local ancestry & Epigenetic Methylation for differential Regions

<p align="center">
  <img src="./imgs/logo.png" alt="LEMUR" width="400"/>
</p>

Stream-merge [modkit](https://github.com/nanoporetech/modkit) **bedMethyl** files from phased haplotypes into a single TSV matrix, with optional local beta-binomial imputation. Built for cohorts where each sample has `*_hp1.bedmethyl` and `*_hp2.bedmethyl` pairs.

Three command-line tools:

| Tool | Purpose |
|------|---------|
| `merge_bedmethyl` | Merge haplotype bedMethyl pairs into one TSV |
| `impute_methylation` | Impute missing methylation on an already-merged TSV |
| `evaluate` | Hold-out benchmark for imputation (per-sample or per-haplotype) |

## Output format (merge)

Tab-separated values (TSV). The first two columns are shared across all samples:

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `pos` | Start position (bedMethyl start, 0-based) |

For each sample label `{id}` given on the command line, columns are appended **in this order**:

#### Haplotype mode (default)

Six columns per sample:

| Column | Source (bedMethyl) | Description |
|--------|-------------------|-------------|
| `{id}.hap1_counts` | `N_modified` (haplotype 1) | Methylated read count |
| `{id}.hap2_counts` | `N_modified` (haplotype 2) | Methylated read count |
| `{id}.hap1_cov` | valid coverage (haplotype 1) | Coverage at the locus |
| `{id}.hap2_cov` | valid coverage (haplotype 2) | Coverage at the locus |
| `{id}.hap1_percentage` | `percent_modified` (haplotype 1) | Methylation 0–100 |
| `{id}.hap2_percentage` | `percent_modified` (haplotype 2) | Methylation 0–100 |

#### Sample mode (`--sample`)

Three columns per sample (hp1 + hp2 aggregated where coverage passes `-c`):

| Column | Description |
|--------|-------------|
| `{id}.counts` | Sum of methylated reads across haplotypes |
| `{id}.cov` | Sum of coverage across haplotypes |
| `{id}.percentage` | `100 × counts / cov` |

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

## Output format (imputed)

After imputation (`--impute` on `merge_bedmethyl`, or `impute_methylation`), the merge columns are replaced by imputed columns. Two output modes are available:

#### Fraction mode (default)

Imputes methylation fraction with a local beta-binomial model.

**Haplotype mode (`--hap`)**

| Column | Description |
|--------|-------------|
| `{id}.hap1_frac_imputed` | Imputed methylation fraction for haplotype 1 (0–1) |
| `{id}.hap2_frac_imputed` | Imputed methylation fraction for haplotype 2 (0–1) |

**Sample mode (default)**

| Column | Description |
|--------|-------------|
| `{id}.frac_imputed` | Imputed methylation fraction for the sample (0–1) |

#### Counts/coverage mode (`--counts-cov`)

Imputes methylated read count and coverage. Fraction is estimated with the same beta-binomial model; coverage is the mean of valid neighbors in the window; counts = round(fraction × coverage).

**Haplotype mode (`--hap`)**

| Column | Description |
|--------|-------------|
| `{id}.hap1_counts_imputed` | Imputed methylated read count (haplotype 1) |
| `{id}.hap1_cov_imputed` | Imputed coverage (haplotype 1) |
| `{id}.hap2_counts_imputed` | Imputed methylated read count (haplotype 2) |
| `{id}.hap2_cov_imputed` | Imputed coverage (haplotype 2) |

**Sample mode**

| Column | Description |
|--------|-------------|
| `{id}.counts_imputed` | Imputed methylated read count |
| `{id}.cov_imputed` | Imputed coverage |

Example header and rows (from `tests/expected/tiny_hap1_imputed.tsv`):

```
chr	pos	S1.hap1_frac_imputed	S1.hap2_frac_imputed
chr1	100	0.4	1
chr1	120	0.25	0
chr1	140	0.5	0.4
```

### Imputation logic

- **Beta-binomial model** with configurable prior (`-a`, `-b`; default uniform 1, 1).
- **Local window** (`-w`, default 200 bp, same chromosome): uses neighboring sites with valid counts/coverage.
- **Minimum neighbors** (`-n`, default 5): imputation runs only when enough valid neighbors exist in the window.
- **Fallback**: if neighbors are insufficient but the site has coverage, writes observed values (fraction, or counts/cov in `--counts-cov` mode); otherwise `.`.
- **Fraction formatting**: `0` and `1` without decimals; other values up to 4 decimal places with trailing zeros stripped (e.g. `0.3333`, `0.425`).

Imputation streams line by line; memory scales with window size × number of haplotype columns, not file size. Use `-j N` to process samples in parallel (OpenMP); each sample’s hap1/hap2 windows are independent.

## Requirements

- C++17 compiler (g++ ≥ 7, clang ≥ 5)
- CMake ≥ 3.14
- OpenMP (optional; enables parallel imputation with `-j`)

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
# or
make
```

Binaries: `build/merge_bedmethyl`, `build/impute_methylation`, `build/evaluate`

```bash
cmake --install build   # optional, installs to CMAKE_INSTALL_PREFIX/bin
```

## Usage

### `merge_bedmethyl`

```bash
merge_bedmethyl [-c N] [-s M] [--sample] [--impute] [-w BP] [-a A] [-b B] [-n N] [-j N] \
  <output.tsv> <label1> <hp1> <hp2> [<label2> <hp3> <hp4> ...]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-c`, `--min-cov N` | Include haplotype fields only if valid coverage (column 9) **>** N | `3` |
| `-s`, `--min-samples M` | Minimum samples with data per row | `N-1` (N = number of pairs) |
| `--sample` | Aggregate hp1+hp2 into `{id}.counts`, `.cov`, `.percentage` | off (haplotype columns) |
| `--impute` | After merge, run beta-binomial imputation (see [Output format (imputed)](#output-format-imputed)) | off |
| `--counts-cov` | With `--impute`: impute counts and coverage instead of fraction | off |
| `-w` | Imputation genomic window (bp, same chromosome) | `200` |
| `-a`, `-b` | Beta-binomial prior α and β | `1`, `1` |
| `-n` | Minimum valid neighbors in window to impute | `5` |
| `-j` | Parallel imputation by sample (`0` = all cores) | `1` |
| `-h`, `--help` | Show help | |

With `--impute`, merge writes a temporary TSV internally, imputes all sample haplotypes, and writes the final imputed matrix to `<output.tsv>`.

#### Example

```bash
./build/merge_bedmethyl merged.tsv \
  CHI01 results/modkit/CHI01.bed/CHI01_hp1.bedmethyl results/modkit/CHI01.bed/CHI01_hp2.bedmethyl \
  CHI02 results/modkit/CHI02.bed/CHI02_hp1.bedmethyl results/modkit/CHI02.bed/CHI02_hp2.bedmethyl
```

Require all samples at each site:

```bash
./build/merge_bedmethyl -s 2 -c 3 merged.tsv CHI01 ... CHI02 ...
```

Merge and impute in one step:

```bash
./build/merge_bedmethyl --impute -w 200 -n 5 -j 4 imputed.tsv \
  CHI01 results/modkit/CHI01.bed/CHI01_hp1.bedmethyl results/modkit/CHI01.bed/CHI01_hp2.bedmethyl \
  CHI02 results/modkit/CHI02.bed/CHI02_hp1.bedmethyl results/modkit/CHI02.bed/CHI02_hp2.bedmethyl
```

### `impute_methylation`

Local **beta-binomial imputation** on an already-merged TSV.

```bash
# Sample mode (default): one column set per sample
impute_methylation [-w 200] [-a 1] [-b 1] [-n 5] [-j N] merged.tsv imputed.tsv

# Haplotype mode: phased hp1/hp2 columns per sample
impute_methylation --hap [-w 200] [-a 1] [-b 1] [-n 5] [-j N] merged.tsv imputed.tsv

# Counts/coverage mode: impute methylated counts and coverage
impute_methylation --counts-cov merged.tsv imputed_counts_cov.tsv
```

| Option | Description | Default |
|--------|-------------|---------|
| `--hap` | Input has `{id}.hap1_counts` / `{id}.hap2_counts` columns | off (`{id}.counts`) |
| `--counts-cov` | Impute counts and coverage instead of fraction | off |
| `-w`, `-a`, `-b`, `-n`, `-j` | Same as merge `--impute` options | `200`, `1`, `1`, `5`, `1` |

Fraction mode (default) drops `{id}.counts`, `.cov`, `.percentage`; writes `{id}.frac_imputed`.
Counts/cov mode (`--counts-cov`) writes `{id}.counts_imputed` and `{id}.cov_imputed`.
Haplotype mode uses the same patterns with `{id}.hap{1,2}_*` prefixes.
Other columns (e.g. `chr`, `pos`) are preserved.

### `evaluate`

Hold-out benchmark with mask-and-impute scoring.

```bash
# Single target, sample mode (default)
evaluate -c CHI08A.counts [-chr chr1] [-m 0.2] [-s 42] [-w 200] [-a 1] [-b 1] [-n 5] merged.tsv

# Single target, haplotype mode
evaluate --hap -c CHI08A.hap1_counts [-chr chr1] ... merged.tsv

# Cohort mode: evaluate all sample/haplotype columns in parallel
evaluate -o cohort.eval.tsv [-chr chr1] [-m 0.2] [-s 42] [-w 200] [-n 5] [-j N] merged.tsv
evaluate --hap -o cohort.eval.tsv ... merged.tsv
```

| Option | Description | Default |
|--------|-------------|---------|
| `--hap` | Input has `{id}.hap1_counts` / `{id}.hap2_counts` columns | off (`{id}.counts`) |
| `-c` | Counts column for single-target mode | — |
| `-o` | Cohort summary TSV (required without `-c`) | — |
| `-chr` | Chromosome to mask, impute, and score | `chr1` |
| `-m` | Fraction of valid sites to mask (hold-out) | `0.2` |
| `-s` | RNG seed for reproducible mask | `42` |
| `-w`, `-a`, `-b`, `-n`, `-j` | Same as `impute_methylation` | `200`, `1`, `1`, `5`, `1` |

Workflow:

1. Masks a reproducible fraction of valid sites in the chosen counts/cov/percentage columns on `-chr`.
2. Writes `<input>.masked.tsv` (single-target) or per-target sidecars (cohort).
3. Imputes the masked column(s) → `<input>.imputed.eval.tsv` or cohort output.
4. Prints MSE and Pearson correlation on masked sites (stderr).

## Slurm

```bash
export SAMPLES="CHI01 CHI02 CHI03"
export MODKIT_DIR=/path/to/results/modkit
export OUTPUT=merged.tsv
sbatch scripts/run_slurm.sh
```

Optional: `MIN_COV` (default `3`), `MIN_SAMPLES` (default: merge tool’s `N-1` rule).

## Tests

```bash
make test
# or
cd build && ctest --output-on-failure
```

CTest runs:

1. **merge_two_samples** — merges fixture bedMethyl pairs and checks row count.
2. **merge_output_format** — compares output to `tests/expected/merge_two_samples.tsv` (see [Output format (merge)](#output-format-merge)).
3. **impute_stream_tiny** / **impute_output_columns** — imputation on `tests/data/tiny.tsv` vs `tests/expected/tiny_hap1_imputed.tsv`.
4. **evaluate_tiny** — hold-out evaluation metrics on the tiny fixture.

### Continuous integration

On push and pull requests to `main` and `dev`, [`.github/workflows/ci.yml`](.github/workflows/ci.yml) builds with CMake and runs `ctest` on Ubuntu.

## Project layout

```
LEMUR/
├── CMakeLists.txt
├── Makefile              # wrapper around CMake
├── imgs/logo.png
├── include/merge_bedmethyl/
├── include/impute_methylation/
├── src/
│   ├── main.cpp          # merge_bedmethyl entry
│   ├── impute_main.cpp   # impute_methylation entry
│   ├── evaluate_main.cpp # evaluate entry
│   └── impute/           # beta-binomial imputation library
├── tests/data/           # small bedMethyl and TSV fixtures
├── tests/expected/       # golden TSV for output tests
├── .github/workflows/    # GitHub Actions CI
└── scripts/run_slurm.sh
```

## Input assumptions

- bedMethyl files are **sorted** by chromosome and start (modkit default).
- modkit bedMethyl column indices used (0-based): start = 1, valid coverage = 9, `percent_modified` = 10, `N_modified` = 11.

## Author

Gabriel Cabas

## License

MIT — see [LICENSE](LICENSE).
