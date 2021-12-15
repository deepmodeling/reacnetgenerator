FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:f9a3bf71ad91d9f0424c3a045de4b7e8eed39b0c5d32ab5aceb90a8671c4b6b3
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
