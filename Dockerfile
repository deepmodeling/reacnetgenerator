FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:f10c188f5c22ee4ad0cb34f13f51efee894736eb19effeaf8b1e7f3ed0d82849
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
