FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:e74e40b67c19d5d084840222d32a57fe711344c2339c102f623a8adea3670cb2
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
