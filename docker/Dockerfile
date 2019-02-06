FROM continuumio/miniconda3
RUN conda install openbabel rdkit -c openbabel -c conda-forge \
    && pip install reacnetgenerator \
    && conda clean --yes --all
