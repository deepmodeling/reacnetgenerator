FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:34c888f89dbfe8615560ff20d42e1c88929be0d73dfe8b9b3006bcfdbffe087f
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
