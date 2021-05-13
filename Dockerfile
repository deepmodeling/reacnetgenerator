FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:09dd32422c3509c9bb086e4b7c6d81871e4e9f041f93a4e473ac4b83518bc71a
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
