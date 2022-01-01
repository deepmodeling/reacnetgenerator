FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:065040ba5773272a1ebf5b55b9ac60552ff2593290c55f73003c2d952854bb96
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
