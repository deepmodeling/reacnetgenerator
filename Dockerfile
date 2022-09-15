FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:279ea354df30a7b1dca1882734079a17730bc0b5b22aa4bd1bab45097ceb51b8
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
