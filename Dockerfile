FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:7f54a6900840e78076ed6196e67ef714547ce13600028c7189635df76c9abaaf
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
