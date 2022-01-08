FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:46cc7ce7e92cc9137c1a896a6f0e748772777e45a35290f5965586ab106a72ec
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
