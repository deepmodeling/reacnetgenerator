FROM continuumio/miniconda3
RUN conda install openbabel rdkit hmmlearn -c openbabel -c conda-forge \
    && pip install reacnetgenerator \
    && conda clean --yes --all \
    && rm -rf /root/.cache
