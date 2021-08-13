FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:29b82bfef93e1ec3b2fe8941e1d64d998dfce6dbe39e55a62d75e5148a5f4e95
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
