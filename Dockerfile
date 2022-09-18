FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:d33fca424b42b1060d53827be7d3f4fece64d1207bfc3b9f55bea60b86ce0a27
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
