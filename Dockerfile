FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:489c5efd3ed8080fc62e84b643d700d5c2cdbac7a17269100970011bef100db0
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
