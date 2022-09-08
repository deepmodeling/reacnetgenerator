FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:39cbfe8564bd7ffeb7507d151b0719ce936bd563c0415b713e264dbc1d5cfae8
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
