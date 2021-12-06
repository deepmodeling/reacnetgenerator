FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:52fae0354cf02750a4a7e4ad30fe7456816ea670cd84759fbe2b256b2b68cab9
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
