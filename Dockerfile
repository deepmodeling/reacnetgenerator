FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:1c30eabd74affd6648c4672a05dc2175666367f96047aea15afe3aff012f0cd8
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
