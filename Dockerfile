FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:fcbef84b0a47df665de4a1cd3e780b2d0e8d6681ac47a46d800e1fed64a01267
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
