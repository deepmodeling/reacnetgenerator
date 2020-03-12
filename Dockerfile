FROM continuumio/miniconda3:latest@sha256:456e3196bf3ffb13fee7c9216db4b18b5e6f4d37090b31df3e0309926e98cfe2
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
