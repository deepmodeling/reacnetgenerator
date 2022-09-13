FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:a37ab93d8dd28af10e9a90d4c88d9a9184363b05fe24a4033031140dd603941c
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
