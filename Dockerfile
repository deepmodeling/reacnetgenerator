FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:f02c7b9bf6d42191a766b2bc8e2d35755dd89bcea92f41de3441bef47e2be43d
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
