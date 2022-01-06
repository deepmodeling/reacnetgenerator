FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:12161c633946d268895e368fe3d44ac91187365c1e2a3113a983c771ceacdb1c
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
