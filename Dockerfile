FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:126c1e73f3c54996145dcae89ea4297a47f001a23e5b3998a7f40161b06edddd
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
