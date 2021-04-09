FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:ed90ee065fb9cf8e10be0eb9faf06cd3ff8cc4793b29fa0b217b5ce98c27c1de
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
