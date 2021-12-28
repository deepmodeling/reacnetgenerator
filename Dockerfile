FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:7e5f84765ac074a8d3a95c44095637473a96b0f99fc01dba983406a10643fc01
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
