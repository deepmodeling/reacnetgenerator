FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:e72480764d31706c41cf542be82d7abcaff18e640f97c47aebdff6a9c70be04f
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
