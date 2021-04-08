FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:88d4812becc9bcce506fddc113575d284075c42fd6f63194cf5f665f12a951d0
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
