# COPLA2

COPLA2 classifies plasmids into Plasmid Taxonomic Units (PTUs). It computes the average
nucleotide identity (ANI) between each query plasmid and the plasmids of a reference
network built from NCBI RefSeq, inserts the query into the network, refines its partition
and reports the PTU of the query with a score.

COPLA2 can be used through the web server at <https://castillo.dicom.unican.es/copla2/>
or run locally with the scripts of this repository, as described below.

## Requirements

- Linux (x86-64) or macOS (Apple Silicon), from a terminal.
- Conda, Mamba or micromamba, to install the dependencies.
- About 48 GB of disk space for the RefSeq 224 reference database, and about 1 GB of
  memory per plasmid classified in parallel.

## Installation

1. Get the code:

   ```bash
   git clone https://github.com/santirdnd/COPLA2
   cd COPLA2
   ```

2. Create and activate the Conda environment (Python 3.10, BLAST+ 2.16.0, BinDash 2.3,
   graph-tool 2.68, NumPy 1.26):

   ```bash
   conda env create -f copla2.environment.yml -n copla2
   conda activate copla2
   ```

   To recreate exactly the tested environments, use the explicit package lists instead:

   ```bash
   conda create -n copla2 --file copla2.linux-64.spec-file.txt
   ```

   (or `copla2.osx-arm64.spec-file.txt` on macOS with Apple Silicon).

3. Download and prepare the RefSeq 224 reference database:

   ```bash
   bin/copla2_setup_db.py databases/Copla_RS224 --download -n 8
   ```

   The script downloads into `databases/Copla_RS224/` the files listed in its
   `manifest.json` (about 2 GB), checks their MD5 checksums, and builds the Copla database
   `CoplaDB.json` (about 46 GB) from the sequences of the curated plasmids. Files already
   in the directory are not downloaded again, unless their MD5 checksum does not match.
   See [Reference databases](#reference-databases).

4. Check the installation with the test of [Testing](#testing).

## Usage

The commands below classify the plasmids of `plasmids.fna` with 8 processors, writing the
results to `results/`. They are run from the repository directory.

1. Validate the input and split it into one FASTA file per plasmid:

   ```bash
   bin/copla2_prepare_queries.py plasmids.fna work/fna \
       --name-mode prefix --ext fna --gzip-output \
       --max-seq-len 1500000
   ```

   For a file with several plasmids, or with the contigs of a draft plasmid, each header
   must have the form `>Plasmid_1|Contig_1`: the plasmid identifier, a vertical bar, and a
   contig identifier unique within that plasmid (see `test/plasmids.fna`). Review the
   `Notes` column of `work/fna/plasmids.manifest.tsv` before continuing.

2. Build the *psani* database of the query plasmids:

   ```bash
   bin/psani.py makedb work/fna/plasmids.fofn work/db \
       -d plasmids --seq-id-source filename -n 8
   ```

3. Compute the ANI between the query plasmids and the reference plasmids:

   ```bash
   mkdir -p results
   bin/psani.py dist work/db/plasmids.json databases/Copla_RS224/CoplaDB.json \
       -o results/plasmids.ani.tsv -n 8
   ```

4. Assign each query plasmid to a PTU:

   ```bash
   bin/copla2.py work/db/plasmids.json databases/Copla_RS224 \
       results/plasmids.ani.tsv results -n 8
   ```

The predictions are in `results/ptu_predictions.tsv`, with one row per query plasmid (PTU,
score and notes). For each query, `results/` also contains `<query>.summary.txt`,
`<query>.related_plasmids.tsv` (the reference plasmids used for the prediction) and, when
the partition is refined, `<query>.sHSBM.tsv`.

`copla2_prepare_queries.py`, `psani.py makedb`, `psani.py dist` and `copla2_setup_db.py`
do not overwrite existing outputs: remove them, or add `--force`. `copla2.py` overwrites
the files of its output directory. Use `-h` with any script for all its options.

## Testing

`test/plasmids.fna` contains four plasmids: ColE1, R388 (split into two contigs, as a
draft plasmid), pXF868 and pXF1.3 (split into two contigs shorter than the 1,000 bp
fragments). Run the commands of [Usage](#usage) with `test/plasmids.fna` as input (a few
minutes with 8 processors) and compare the `Predicted_PTU` column of
`results/ptu_predictions.tsv` with that of `test/expected/<platform>/ptu_predictions.tsv`,
where `<platform>` is `linux-64` or `osx-arm64`:

| Plasmid | Expected PTU |
|---|---|
| NC_001371.1 (ColE1) | PTU-E1 |
| NC_028464.1 (R388) | PTU-W |
| NC_002092.1 (pXF868) | none (no ANI link to the reference plasmids) |
| NC_002489.3 (pXF1.3) | none (no fragments to compare with the reference plasmids) |

With the default seed (`--seed 41`), the scores also match on the same platform and with
the same software versions. The ANI values can also be compared with
`test/expected/<platform>/plasmids.ani.tsv`, sorting both files, since the order of their
rows changes between runs.

## Reference databases

Each reference database is a directory of `databases/` with a `manifest.json`, the only
file of the directory kept in the repository. The manifest lists the files of the
database, with their MD5 checksum and download URL, and the parameters used to build the
Copla database and the reference network.

Each version of a database has its own directory (e.g. `databases/Copla_RS224_v1.0/`), so
several versions can be installed side by side, and `databases/Copla_RS224` is a link to
the latest one, used in the commands of this README. When a new version is released,
prepare it with `copla2_setup_db.py databases/Copla_RS224 --download`; the directory of
the previous version can be kept, to reproduce earlier results, or deleted.

For RefSeq 224 version 1.0 (`databases/Copla_RS224_v1.0/`):

| File | Content |
|---|---|
| `README.md` | Description of the database |
| `RS224_plasmid_sequences.zip` | Sequences of the 89,393 RefSeq 224 plasmids (`fna/<ID>.fna.gz`) |
| `RS224_plasmid_curated.lst` | Identifiers of the 64,030 plasmids of the reference network |
| `RS224_plasmid_network.pickle` | Reference network, with its partition and the PTU of each plasmid |
| `RS224_plasmid_metadata.tsv.xz` | Metadata of every plasmid, including the curation and its PTU (optional) |
| `RS224_plasmid_network.ani.tsv.xz` | ANI between the plasmids of the network, in the format of `psani.py dist` (optional) |

The optional files are not needed by the COPLA2 scripts: `copla2_setup_db.py` verifies
them when present, but downloads them only with `--download-optional` (instead of
`--download`).

`copla2_setup_db.py` builds from them `CoplaDB.json`, with the fragment databases (`fwa/`)
and k-mer sketches (`kmer/`) of the curated plasmids. The paths in `CoplaDB.json` are
relative to the database directory, so the database can be moved without being prepared
again (`--absolute-paths` writes absolute paths instead). The per-plasmid FASTA files are
deleted once the database is built (`--keep-fna` keeps them). If the files have been
downloaded beforehand, place them in the database directory and run the script without
`--download`.

The database files are archived in Zenodo:

| Database | Directory | DOI |
|---|---|---|
| RefSeq 224, version 1.0 | `databases/Copla_RS224_v1.0/` | [10.5281/zenodo.23989163](https://doi.org/10.5281/zenodo.23989163) |
| RefSeq 200, version 1.0 | `databases/Copla_RS200_v1.0/` | [10.5281/zenodo.23990348](https://doi.org/10.5281/zenodo.23990348) |

## Scripts

| Script | Purpose |
|---|---|
| `bin/copla2_setup_db.py` | Downloads and prepares a reference database described by its `manifest.json` |
| `bin/copla2_prepare_queries.py` | Validates the input FASTA and splits it into one file per plasmid |
| `bin/psani.py` | Builds fragment databases (`makedb`) and computes ANI between them (`dist`) |
| `bin/copla2.py` | Inserts the query plasmids into the reference network and assigns their PTU |

*psani* re-implements the BLASTn-based two-way ANI of the Enveomics collection, runs the
comparisons in parallel and skips the plasmid pairs with no k-mer similarity, estimated
with BinDash.

## Citation

COPLA:\
Redondo-Salvo S, Bartomeus-Peñalver R, Vielva L, et al. (2021) COPLA, a taxonomic
classifier of plasmids. BMC Bioinformatics 22:390.
https://doi.org/10.1186/s12859-021-04299-x

PTUs:\
Redondo-Salvo S, Fernández-López R, Ruiz R, et al. (2020) Pathways for horizontal gene
transfer in bacteria revealed by a global map of their plasmids. Nat Commun 11:3602.
https://doi.org/10.1038/s41467-020-17278-2

## License

COPLA2 is distributed under the GNU General Public License v3.0 (see `LICENSE`).
`psani.py` is a re-implementation of `ani.rb` from the Enveomics collection
(<https://github.com/lmrodriguezr/enveomics>, Artistic License 2.0).
