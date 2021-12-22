FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:89bde6bf4fd0ee6fdb5484a54f6c19daf9cc3a74c6b81df46faa2f2e8f89896c
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
