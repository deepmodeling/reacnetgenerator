FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:2d1223fbf77bfdc97895e16c5f9355c962820bc3b0311e5f1ceeb159c75ba1be
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
