FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:e2ec7a99132452dfcf32c5ad39b4c73b2fcc4aba7f6000f5dc3c77444dd25566
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
