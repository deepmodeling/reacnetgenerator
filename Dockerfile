FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:945e3cd0649ef6af5304a5c3fbeb082b6e16c980a55a46d16b1cc66c571d4483
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
