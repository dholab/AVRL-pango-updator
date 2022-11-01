FROM condaforge/mambaforge

# setting the shell to bash, not sh
SHELL ["/bin/bash", "-c"]

# set channels for conda
RUN conda config --add channels defaults \
&& conda config --add channels bioconda \
&& conda config --add channels r

# install python
RUN mamba install -y python && conda clean -afy

# install analysis tools
RUN mamba install -y iqtree && conda clean -afy
RUN mamba install -y -c bioconda pangolin && conda clean -afy
RUN mamba install -y -c bioconda bbmap && conda clean -afy
RUN mamba install -y -c conda-forge curl && conda clean -afy
RUN mamba install -y -c conda-forge zstd && conda clean -afy
RUN mamba install -y -c conda-forge r-base && conda clean -afy
RUN mamba install -y -c bioconda minimap2 && conda clean -afy
RUN mamba install -y -c bioconda snpeff && conda clean -afy
