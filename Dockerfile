FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:77dac2d889ecec601693c4a0db2b01fb974d835c671955ac250e2bb4f9aa74f0
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
