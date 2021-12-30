FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:866b5ac300fd109055f9a989f12901479017d3aa657f4ce3e8bca3ac472445e9
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
