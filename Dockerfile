FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:b2c93c8be04c78989ec753a0fffe8cfa839a14a4e54144f48a17b5a916fb9ae7
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
