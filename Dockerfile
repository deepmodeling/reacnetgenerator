FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:0c2628283db28126a6614fccfa695734d2ed8b853dcfb69a87bf502969050786
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
