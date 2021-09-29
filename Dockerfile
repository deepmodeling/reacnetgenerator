FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:dcae66b5c116f27be19e1950828d99510617d7d7ac372f725177a9072f0a5761
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
