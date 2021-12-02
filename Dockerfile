FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:6a636779a7c08149796e79e10d99485b2f051fc3263a1c549dd63d18bbbb8059
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
