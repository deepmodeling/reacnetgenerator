FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:d6fde81b93b9e2d10d2956d52ebbdc1893ae5a4e99be39c1c80009cbd0b64125
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
