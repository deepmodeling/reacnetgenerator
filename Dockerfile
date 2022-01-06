FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:d90f939b31e9509ec24072fbf75fb0c48f5b6fdd244b705f6df53373411a1728
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
