FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:24d17ab5340a92d5c958cb4e4bb9de1f03800fa1b9e938ca9ddfdaf9e8b06e90
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
