FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:dad75e757dba79368427589011ff850652d85d5c53c1828b1cb13e5cbc845fdd
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
