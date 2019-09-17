FROM continuumio/miniconda3:latest@sha256:54eb3dd4003f11f6a651b55fc2074a0ed6d9eeaa642f1c4c9a7cf8b148a30ceb
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
