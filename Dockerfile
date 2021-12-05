FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:9713bbd55befc89808ff3749ac452931698917e02fecba62c9bdad3fc4b79487
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
