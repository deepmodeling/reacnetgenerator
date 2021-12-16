FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:24113fa6bb367f1524b4ca3ad807ae04205039b76fea4e9294e7d280ce8a0f59
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
