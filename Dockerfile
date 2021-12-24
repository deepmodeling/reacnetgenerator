FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:13bcd40a6931f7369bb78f39bdacfcaf5651c1e2d715d8695ec5274a3b433f5e
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
