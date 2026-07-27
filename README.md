# SpliceImpactR Studio

SpliceImpactR Studio is a browser-based Shiny application for connecting
alternative-splicing changes to transcript, protein, domain, and
protein-protein interaction consequences.

## Use the web app

**[Launch SpliceImpactR Studio](https://fiszbein-lab.shinyapps.io/SpliceImpactR_Studio/)**

The hosted app provides a guided demonstration as well as workflows for your
own annotations, splicing evidence, protein features, and interaction data.

### Quick demonstration

1. Open the app and click **Load Demo Workspace**.
2. Confirm the comparison. The demonstration defaults to reference
   **control** and case **case**.
3. Click **Run Full Pipeline**. The app loads the selected bundled PPI network
   and runs differential inclusion, transcript mapping, sequence/frame
   comparison, domain analysis, and PPI-switch analysis.
4. Review the **DI**, **Consequences**, **Explore**, **Enrichment**, and
   **Exports** tabs.

The demonstration uses a small GENCODE 45 subset and packaged example data. It
is intended for learning and validation, not genome-wide analysis.

## Analyze your own data

The controls on the left follow the workflow below.

### 1. Load a reference

Choose one of the available annotation sources and click **Load Annotations**:

- **Bundled human GENCODE 45 reference — recommended:** uses the full,
  preprocessed GRCh38.p14 reference supplied by the hosted app. Nothing needs
  to be downloaded or uploaded.
- **Bundled test set:** loads the small reference used by the demonstration.
- **Load GENCODE 45 files from this computer:** upload the official annotation
  GTF, transcript FASTA, and protein FASTA files. The app provides direct
  GENCODE download links for obtaining these files on your computer first.

Browser-uploaded reference files are used only in an isolated app session and
are removed when that session ends.

### 2. Load splicing evidence

Open **Event Inputs** and select one of the following:

- a self-contained rMATS project ZIP containing its sample manifest;
- a normalized raw-event table with sample-level evidence; or
- a post-differential-inclusion event table.

After loading the data, confirm which condition is the reference and which is
the case. When conditions named `control` and `case` are present, the app
selects them automatically.

### 3. Load protein features

Open **Protein Features** and either:

- query matching features from Ensembl/BioMart; GENCODE 45 corresponds to
  Ensembl 111;
- upload a normalized SpliceImpactR feature table; or
- upload a table containing manual protein coordinates.

The packaged feature subset is compatible only with the demonstration
reference. A full-reference analysis requires queried or uploaded features.

### 4. Select interaction data and run

Under **PPI + Background**, keep the bundled SpliceImpactR interaction network
or upload a normalized PPI table. **Run Full Pipeline** loads the selected
network on demand and performs the complete analysis. Individual stages can
also be run from **Analysis Controls**.

## Results and downloads

The app provides interactive previews and plots for differential inclusion,
matched transcript pairs, nucleotide and protein consequences, domain changes,
PPI changes, transcript exploration, and gene-set enrichment. The **Exports**
tab downloads result tables as CSV files and also provides an activity log and
provenance record.

## Project

SpliceImpactR Studio is an interactive interface for the
[fiszbein-lab/SpliceImpactR](https://github.com/fiszbein-lab/SpliceImpactR) R
package. This lightweight repository contains the application code; it does
not publish user data, reference bundles, cached downloads, or analysis
outputs.

For local development, restore the locked dependencies and launch from the
repository root:

```r
renv::restore()
shiny::runApp(".")
```

Run the automated checks with:

```sh
Rscript tests/testthat.R
Rscript tests/smoke_demo.R
```
