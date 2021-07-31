FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:6fc21664e3170b7eb39a8a8167e5283390fadbf2c04a6637e94a61a97d622447
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
