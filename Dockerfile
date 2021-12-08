FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:d2684b3f36680e171e70b3f85ad9ce7632548dc50095052430d9fc3fe4c46e65
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
