FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:d8c0c06f48194cfb97f68c6ef96ea15194c0ce31077bd7c46c4210d2c39291d2
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
