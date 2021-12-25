FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:35bd80fb95b5e0440da8242496ae155ba4a0844ae5b954e0c77f2b5390828853
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
