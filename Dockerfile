FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:ee73c81544ef9c62094986dd89bb3dc2d11cb1fdb262ee6329b10733d8e82bc2
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
