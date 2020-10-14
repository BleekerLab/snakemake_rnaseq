FROM continuumio/miniconda:4.7.12

WORKDIR /home/snakemake/

COPY environment.yaml ./


# mamba is a faster C++ re-implemenbtation of conda
RUN conda install -c conda-forge mamba --yes \
  && mamba env create -f environment.yaml \
  && conda clean --all

RUN echo "source activate rnaseq" > ~/.bashrc
ENV PATH /opt/conda/envs/rnaseq/bin:$PATH

ENTRYPOINT ["snakemake"]
