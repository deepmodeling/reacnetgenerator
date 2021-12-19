FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:1de0eba99cbb2b5f169e74040ecf18a99e704260260630f16f5a792dc1e3b070
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
