FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:6719bd6a073ed9dd856a717a6c251dfb7de6c5f5806d25dbe10e8ea7f8385b1b
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
