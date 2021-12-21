FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:95c380d8a112b93175d531d4dfd7ef71280e45f0e9c14feeaa650018747c3892
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
