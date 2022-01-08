FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:4d8759406076d9af6c8e016fb1c2b2ce53fbd607a469aff364c17f4e921a858f
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
