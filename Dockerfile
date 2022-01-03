FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:f7769697d8c5305e0f7af0809a1615f1852b64a3c48cd97b81de00437c3c3e2d
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
