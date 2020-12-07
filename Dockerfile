FROM continuumio/miniconda3:latest@sha256:7838d0ce65783b0d944c19d193e2e6232196bada9e5f3762dc7a9f07dc271179
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
