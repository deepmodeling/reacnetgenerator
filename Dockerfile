FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:f8ea772d3343bdec2b23f3280a5f358b6b5759d5f556290a22bf9375f18c08ce
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
