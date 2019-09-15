FROM continuumio/miniconda3:latest
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
