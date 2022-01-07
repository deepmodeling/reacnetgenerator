FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:d23497ca9ce7e4a52b8f6b54f59452ff28b51173f6f95bbbb0cd2435da154124
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
