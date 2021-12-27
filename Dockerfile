FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:e20efdb1aa6676e18fa2ce0adcecf9cf1f106e325c3de73a84dbd810554f2694
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
