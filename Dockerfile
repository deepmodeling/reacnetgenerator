FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:ed9e9bcc0a186c4c23d6b2163528475eb6c1b63b44635c3d0e110bde5c863829
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
