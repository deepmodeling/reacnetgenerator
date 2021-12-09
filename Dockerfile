FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:0c06b7917ff4daba9ce0482e40b8870b591fc190bf113eaf69e601cbd8a40837
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
