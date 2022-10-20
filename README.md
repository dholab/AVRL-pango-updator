# NextFlow Pipeline for Analyzing SARS-CoV-2 Pangolin Lineages

## Overview

This pipeline was designed to update Pangolin lineages for an arbitrary number of SARS-CoV-2 consensus sequences and then extend these findings with some additional analysis. In all, the pipeline goes through three major steps

1. Ensuring that it uses the most updated version of [Pangolin](https://github.com/cov-lineages/pangolin) to reclassify any number of FASTA-formatted consensus sequences. _NOTE:_ This pipeline does not assemble consensus sequences itself. These files must be produced beforehand.
2. Identifying putative prolonged infections from the now-reclassified consensus sequences.
3. Classifying "RBD Mutation Level," a metric for how evolutionarily advanced an infection is, for each sample sequenced since the emergence of lineage BA.2.

## Quick Start

If you already have miniconda3 and NextFlow installed on your system, simply invoke the pipeline with:

```
nextflow run dholab/AVRL-pango-updator \
--update_pango true \
--identify_long_infections true \
--classify_mutation_levels true \
--days_of_infection 240 \
--data_dir "/abolute/path/to/your/data/files/" \
--results results/ \
-profile conda
```

Alternatively, if you prefer to use Docker, invoke the workflow with the Docker profile (_NOTE: This arm of the pipeline is still under development_), like so:

```
nextflow run dholab/AVRL-pango-updator \
--update_pango true \
--identify_long_infections true \
--classify_mutation_levels true \
--days_of_infection 240 \
--data_dir "/abolute/path/to/your/data/files/" \
--results results/ \
-profile docker
```

If you do not have NextFlow, miniconda3 and/or Docker installed on your system, or you want more information about the pipeline's configuration, proceed to the following steps.

## Detailed Setup Instructions

### Inputs and Outputs

### Pipeline Configuration

To invoke, run:

```
nextflow run dholab/AVRL-pango-updator \
--update_pango true \
--identify_long_infections true \
--classify_mutation_levels true \
--days_of_infection 240 \
--data_dir "/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data" \
--results results/ \
-profile conda
```

Or, if you would prefer to clone the workflow first:

```
git clone https://github.com/dholab/AVRL-pango-updator.git .

nextflow run main.nf \
--update_pango true \
--identify_long_infections true \
--classify_mutation_levels true \
--days_of_infection 240 \
--data_dir "/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data" \
--results results/ \
-profile conda
```
