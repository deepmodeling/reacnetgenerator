FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:5f3f6ec135d2bc7d5740a9a23ea06c62427c7aaf4034207cdc7e1d55a5129788
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
