FROM condaforge/mambaforge

# setting the shell to bash, not sh
SHELL ["/bin/bash", "-c"]

# set channels for conda
RUN conda config --add channels defaults \
&& conda config --add channels bioconda \
&& conda config --add channels r

# install python and related utilities
RUN mamba install -y python && conda clean -afy
RUN mamba install -y -c anaconda pip && conda clean -afy
RUN mamba install -y -c conda-forge biopython && conda clean -afy

# install conda-based analysis tools
RUN mamba install -y -c bioconda bbmap && conda clean -afy
RUN mamba install -y -c conda-forge curl && conda clean -afy
RUN mamba install -y -c conda-forge zstd && conda clean -afy
RUN mamba install -y -c conda-forge r-base && conda clean -afy
RUN mamba install -y -c bioconda snakemake && conda clean -afy
RUN mamba install -y -c bioconda minimap2 && conda clean -afy
RUN mamba install -y -c bioconda snpeff && conda clean -afy

# transferring pangolin scripts
ADD pangolin/ pangolin/

# install pangolin and dependencies with pip
RUN cd pangolin && python setup.py install
RUN pip install snakemake
RUN pip install git+https://github.com/cov-lineages/pangolin-data.git
RUN pip install git+https://github.com/cov-lineages/scorpio.git
RUN pip install git+https://github.com/cov-lineages/constellations.git
