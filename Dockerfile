FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:b24e656d8b4150a72aaeaaea58063b3ae5634b289f8547d7934d229589ab6718
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
