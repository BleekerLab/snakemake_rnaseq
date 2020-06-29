FROM continuumio/miniconda:4.7.12

WORKDIR /home/snakemake/

COPY environment.yaml ./

RUN conda env create -f environment.yaml

RUN echo "source activate rnaseq" > ~/.bashrc
ENV PATH /opt/conda/envs/rnaseq/bin:$PATH

ENTRYPOINT ["snakemake"]