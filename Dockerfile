FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:89665dd077939bad9f76df287d1a91cd7a3b66f50570281e4f07731679654e76
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
