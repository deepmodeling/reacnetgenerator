FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:cf449cdbe86c00aee02f89c5abf9755ee4f8f87c9d250421c6e1709625a4d207
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
