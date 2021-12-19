FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:529d4a08c41dd68fc9f2176dc681df91cf1588451dff4443043d52bb1741939c
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
