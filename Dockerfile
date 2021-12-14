FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:18d7137b52a40c1cb66ef6afa683d033ba2b808458f14d55ef72413395e9f8a1
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
