FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:f31450e6ff28f67cf9fcb0d533b0282d5934a8b8ee8ea68a377d47c6c59a3e08
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
