FROM condaforge/mambaforge

# set channels for conda
RUN conda config --add channels defaults \
&& conda config --add channels conda-forge \
&& conda config --add channels bioconda \
&& conda config --add channels r

# install analysis tools
RUN mamba install -y iqtree && conda clean -afy
RUN mamba install -y pangolin && conda clean -afy
RUN mamba install -y bbmap && conda clean -afy
RUN mamba install -y curl && conda clean -afy
RUN mamba install -y zstd && conda clean -afy
RUN mamba install -y r-base && conda clean -afy
RUN mamba install -y minimap2 && conda clean -afy
RUN mamba install -y snpeff && conda clean -afy
