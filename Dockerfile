FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:2f052c9ecbd7cd211868a087d1ef56b6cd106836083a3935866e34c21a78fa3a
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
