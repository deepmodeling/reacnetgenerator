FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:3291382e3f70f02dfc18033a4cbdaac55c1620ea4820bff38addc4cdb3cedd84
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
