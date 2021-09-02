FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:948a203e999a8988f22a6c29eead69ac987bd3eba9f97d11fdaecc9f71ae54c5
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
