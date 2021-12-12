FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:dcc83322ab9ccc84511d385f79500f8fbb83d84b3f8ef8e259179c9bcc79986c
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
