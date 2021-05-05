FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:3cda15fe15a4d966093a409fdec1392c1b3e9fe788bb0d72b2e6f50a54d41b1d
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
