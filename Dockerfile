FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:5fa57fb5216d902299079c9b1a7ee417fdc8b42dc6899ab5949be9c10655b379
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
