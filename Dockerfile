FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:18da5f25134ce43173521546bcb3b62a3b20697ba6ae8f959756e7ccdc775551
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
