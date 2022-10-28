# NextFlow Pipeline for Analyzing SARS-CoV-2 Pangolin Lineages

## Overview

This pipeline was designed to update Pangolin lineages for an arbitrary number of SARS-CoV-2 consensus sequences and then extend these findings with some additional analysis. In all, the pipeline goes through three major steps:

1. Ensuring that it uses the most updated version of [Pangolin](https://github.com/cov-lineages/pangolin) to reclassify any number of FASTA-formatted consensus sequences. _NOTE:_ This pipeline does not assemble consensus sequences itself. These files must be produced beforehand.
2. Identifying putative prolonged infections from the now-reclassified consensus sequences.
3. Classifying "RBD Mutation Level," a metric for how evolutionarily advanced an infection is, for each sample sequenced since the emergence of lineage BA.2 (learn more about this lineage on [Open cov-Spectrum](https://open.cov-spectrum.org/explore/United%20States/AllSamples/AllTimes/variants?pangoLineage=BA.2&) and on [outbreak.info](https://outbreak.info/situation-reports?pango=BA.2)).

For each of these three steps, the pipeline has two modes: An "all sequences mode" and a "runs of interest mode." If no runs of interest are specified in the command line or configuration file (more details below), it will default to finding and processing all FASTA files in the provided file path. If runs of interest are provided, it will only run on those runs' sequences.

## Quick Start

If you already have Docker and NextFlow installed on your system, simply invoke the pipeline with:

```
nextflow run dholab/AVRL-pango-updator -latest \
--data_dir /abolute/path/to/your/data/files/ \
--results /abolute/path/to/your/results/
```

Alternatively, you can also specify which analyses you would like to skip, like so:

```
nextflow run dholab/AVRL-pango-updator -latest \
--update_pango false \
--identify_long_infections false \
--classify_mutation_levels false \
--data_dir /abolute/path/to/your/data/files/ \
--results /abolute/path/to/your/results/
```

If you do not have NextFlow, miniconda3 and/or Docker installed on your system, or you want more information about the pipeline's configuration, proceed to the detailed instructions after the AVRL usage section.

#### A Note on In-House Usage at AVRL

This pipeline was originally developed for use at the [AIDS Vaccine Research Laboratory (AVRL)](https://dholk.primate.wisc.edu/wiki/home/page.view?name=home_index) at the University of Wisconsin-Madison. If you are on staff at AVRL, the pipeline can be invoked using its defaults with a simpler command:

To run the pipeline on all available sequences and perform all analyses, simply run:

```
nextflow run dholab/AVRL-pango-updator -latest
```

To run the pipeline on specific sequencing runs, use this command with the "DHO\_" experiment numbers of interest in a comma-separated list:

```
nextflow run dholab/AVRL-pango-updator -latest --runs_of_interest DHO_12345,DHO_9789
```

## Detailed Setup Instructions

First, make sure you are using a POSIXct-compatible system (e.g. MacOS, Linux, [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install)) with Git installed.

Then, download the pipeline files by running `git clone` in a directory of your choice:

```
git clone https://github.com/dholab/AVRL-pango-updator.git .
```

After the repository has downloaded, you may need to set the pipeline's scripts to executable with `chmod +x bin/*`, though this should already be done.

If you choose to use conda to organize software dependencies with `-profile conda`, we recommend you install the miniconda python distribution, if you haven't already: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

If you choose to use Docker to organize software dependencies, as is the pipeline's default, the Docker engine must be installed. To do so, simply visit the Docker installation page at [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/).

### Nextflow Installation

This pipeline was built with the [NextFlow](https://www.nextflow.io/) pipeline manager (v22.10.0.5826). We recommend you install NextFlow to your system in one of the two following ways:

#### 1) Installation with Conda

1. If you did not install conda in the steps above, go to the following link to download Miniconda: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
2. Install the `mamba` package installation tool in the command line:
   `conda install -y -c conda-forge mamba`
3. Install Nextflow to your base environment:
   `mamba install -c bioconda nextflow `

#### 2) Installation with curl

If you would prefer not to use conda at all, you can also install NextFlow with curl:

1. Run the following line in a directory where you'd like to install NextFlow:
   `curl -fsSL https://get.nextflow.io | bash`
2. Add this directory to your $PATH. If on MacOS, a helpful guide can be viewed [here](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/).

To double check that the installation was successful, type `nextflow -v` into the terminal. If it returns something like `nextflow version 21.04.0.5552`, NextFlow has successfully installed.

With Docker/conda, NextFlow,

### Running the pipeline

The pipeline is invoked from the working directory where you downloaded it with one shell command:

```
nextflow run main.nf \
--update_pango true \
--identify_long_infections true \
--classify_mutation_levels true \
--data_dir /abolute/path/to/your/data/files/ \
--results results/
```

This command will run the pipeline on all available FASTA files ("all sequences mode") and specifies some key information for the pipeline:

- `--update_pango true` tells the pipeline to update Pangolin
- `--identify_long_infections` tells the pipeline to find putative prolonged infections
- `--classify_mutation_levels` tells the pipeline to classify samples according to the RBD mutation level
- `--data_dir` tells the pipeline where to look for input FASTA files. Note that this pipeline is currently written to find FASTA files in sequencing run subdirectories called `gisaid/`, as is the case in the sequencing process at [AVRL](https://dholk.primate.wisc.edu/wiki/home/page.view?name=home_index). This hard-coding will be replaced with something more flexible in future updates.
- `--results` tells the workflow where to place output files. This can be an absolute path to any location on your system, or a relative path within the workflow directory. By default, results will be placed in a subdirectory of the workflow directory called `results`.
- Note also that the pipeline supports using conda or Singularity instead of Docker; do so by adding `-profile conda` or `-profile singularity` to the nextflow run command. Note that the `-profile` flag has only one dash because single-dash flags modify the behavior of NextFlow itself, whereas double-dash flags specify parameters in the configuration file `nextflow.config`. Here, we do not specify a profile, so Nextflow uses Docker as is specified in its standard profile.

If the pipeline crashes for any reason, simply resume it with:

```
nextflow run main.nf \
--update_pango true \
--identify_long_infections true \
--classify_mutation_levels true \
--data_dir "/abolute/path/to/your/data/files/" \
--results results/ \
-resume
```

To run the pipeline on only a few sequencing runs of interest, specify those runs with the `--runs_of_interest` flag, separated by commas:

```
nextflow run main.nf --runs_of_interest run1,run2,run3
```

By only specifying the `--runs_of_interest` flag, the pipeline will perform all analyses by default and will use its default input and results directory settings. These can be found in the file `nextflow.config`.

_NOTE_ that the run names you specified must also be the name of a folder within which the FASTA files can be found.

## Pipeline Compute Resources

This pipeline will run most efficiently when at least eight (logical) cores and 16 gigabytes of RAM, though the more cores you can make available to it, the faster it will run.

## Inputs and Outputs

As an input, specified with the flag `--data_dir`, the pipeline takes a folder to search recursively through for FASTA-formatted sequence files. These sequences must be consensus sequences.

The workflow produces four outputs in the folder specified with `--results`:

1. A report with updated pangolin lineages for all found consensus sequences, along with RBD Mutation Levels for samples sequenced after the emergence of lineage BA.2.
2. A reference fasta for the lineage BA.2.
3. A table of putative long infections from among all detected consensus sequences.
   - These infections come from samples that classify as "old" lineages, i.e. lineages that were first designated in pangolin long before a sequence classified as it in a sequencing run. By default, this amount of time is eight months, but this is subject to change in the future. In such cases, where a lineage appears long after it arose and subsided, it is more likely that the infected individual has sustained a prolonged infection since that lineage was prevalent, and less likely that the old lineage re-appeared despite competition from newer, more fit lineages.

Finally, the workflow also produces a visualization of itself called `lineager-analyzer-visualization.png`. This is analogous to a Directed Acyclic Graph (DAG).

## Pipeline Configuration

This pipeline is configured with parameters in the file `nextflow.config`, many of which have already been mentioned above. In general, file paths and pipeline settings are specified there instead of having them hard-coded into the pipeline script `main.nf`.

The configuration parameters include:

- `data_dir`: Absolute path to the directory where subdirectories for each sequencing run are stored.
- `results`: Where to place results. If the path is relative, the results directory will be created within the workflow directory.
- `refgff`: Path to SARS-CoV-2 annotations. This is not used in the current version of the pipeline but may be implemented in the future.
- `distribute_results`: whether to distribute pango lineage reports for each FASTA (`true` or `false`) into the directories where those FASTA files were found (Default: `false`). The lineage reports will contain the date they were created in the file name, allowing you to compare lineage reports from different dates and different versions of pangolin.
- `identify_long_infections`: whether to identify (`true` or `false`) potential long infections (Default: `true`)
- `days_of_infection`: How many days past lineage designation to consider an infection prolonged (Default: 240)
- `classify_mutation_levels`: whether to classify into BA.2-based RBD mutation levels (`true` or `false`; default is `true`)
  - This process uses an R script to trim down all observed variants to only those in the Spike protein receptor binding domain, i.e. Spike amino acid residues 319–541. It then adds the count of RBD mutations——an RBD mutation "level", a la Cornelius Roemer's method for tracking convergent evolution among SARS-CoV-2 lineages——to the new pango lineage classifications.
- `update_pango`: whether to update pango to the latest version (`true` or `false`; default is `true`)
- `docker_reg`: Docker registry to use. In the past we have used 'dockerreg.chtc.wisc.edu/dabaker3'
