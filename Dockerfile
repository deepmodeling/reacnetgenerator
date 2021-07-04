FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:e169daf90417f7ef88860707b582f807af0ff833a15dc0fa481b2b1793bef898
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
