FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:511ca5d9b3f3f34fb5c3e8bbcdd2b4bf585c6985657a39207a51e260210e845f
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
