FROM continuumio/miniconda3:latest@sha256:1069d3220e8b6c8259d07f3fafbc36aff3d43a02a6d0efac7f07d139b458cf65
COPY . /reacnetgenerator
RUN conda config --add channels conda-forge && \
    conda install conda-build -y && \
    conda build /reacnetgenerator/conda/recipe && \
    conda install reacnetgenerator --use-local -y && \
    conda clean -tipsy && \
    conda build purge-all && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
