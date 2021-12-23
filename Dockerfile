FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:489dfa2e4f5c2b733b9611f79c72fc509121cb0d6734067a3c0927920b6b4c6e
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
