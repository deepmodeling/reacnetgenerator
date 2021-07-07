FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:17dc2eb414c0937dbb8e6a62eba18afea3635aa4ee1040d508f1c572967cc2bc
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
