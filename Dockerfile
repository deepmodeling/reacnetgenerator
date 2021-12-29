FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:1bef34df3a08be4c4eb4c236fcb111e82ea55af979674ec4589e88b31c002c06
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
