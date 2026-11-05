# Demultiplexing Benchmark

## Code Structure

```text
BAM_demultiplex/
|-- config.sh                 # Paths and benchmark settings
|-- environment_setup/        # Optional dependency installation
|-- Data_Simulation/          # Construct pools from saved barcode lists
|-- Reference_Preparation/    # Prepare donor SNP and HLA references
|-- Benchmarking/             # Run selected demultiplexing methods
`-- Evaluation/               # Evaluate saved predictions
```

## Dependencies

Linux/Bash and Python with pandas. Install tools needed for the selected step:

- Pool construction: subset-bam and samtools.
- Reference preparation: bcftools, tabix, samtools, arcasHLA, Cell Ranger and sinto.
- Demultiplexing: HLA-Splitter and its dependencies, CellSNP-lite, Vireo, Popscle, or Souporcell with Apptainer/Singularity.
- Resource logging: GNU `/usr/bin/time`.

Activate the environment manually. `environment_setup/setup_environment.sh` provides an optional Conda setup; install HLA-Splitter, Cell Ranger and the Souporcell image separately.

## Path Settings

Edit `config.sh`: set `PROJECT_DIR`, `REF_GENOME`, `HLADB`, `SOUP_SIF` and `COMMON_SNP_VCF`; review donor IDs, `SAMPLE_NUM`, `CELL_NUMS_STR` and `THREADS`.

```text
PROJECT_DIR/
|-- BAM_file/P1/possorted_genome_bam.bam  # Indexed BAM for each donor
|-- barcode/100cells-1.txt              # Selected cell barcodes
|-- True_sample_labels/100cells-1.csv   # barcode_select and sample columns
|-- Demuxlet_VCFs/                     # Donor VCFs
|-- HLA_Genotypes/                     # Donor HLA table
`-- Extract_file/                      # Pool BAMs and tool outputs
```

`COMMON_SNP_VCF` defaults to `Revison/Benchmark_inputs/common_variants/00-common_all.chr_added.vcf.gz`, outside the code directory. Set `MERGED_DONOR_VCF` and `HLA_LIST` to the reviewed donor references. Genome builds, contig names and donor labels must match the input data.

## Manual Execution

Run commands individually from `BAM_demultiplex/`; inspect each stage before continuing. No command runs the entire workflow automatically.

```bash
# Set paths in config.sh, then load settings in the current shell.
source config.sh

# 1. Construct pools using existing barcode lists.
bash Data_Simulation/01_simulate_pools.sh

# 2. Prepare references; check donor names and HLA calls afterward.
bash Reference_Preparation/02_prepare_references.sh

# 3. Select a pool and method; adjust as needed.
export CELL_NUMS_STR="100cells-1"
export TOOLS="HLA-Splitter"
bash Benchmarking/03_run_benchmark.sh

# 4. Evaluate existing predictions for the selected pool.
python Evaluation/04_evaluate_predictions.py
```

Available `TOOLS`: `HLA-Splitter`, `Vireo_free`, `Vireo_genotype`, `Demuxlet`, `Souporcell` (space-separated). Set `HLA_GRAPH_CORR=1` explicitly if neighborhood refinement is required. Existing nonempty prediction files are skipped unless `RERUN_EXISTING=1`; inspect incomplete outputs before rerunning. Evaluation writes `benchmark_summary.csv` and `.json` in `PROJECT_DIR`, replacing existing summaries.

Only invoke the stages needed for the current analysis; existing pools and references can be reused. There is no workflow dispatcher or automatic full-pipeline command.
