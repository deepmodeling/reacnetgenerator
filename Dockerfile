FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:9f2ccfa16e3c5298896de78e5073382ea14c1aeff46f7aeb6016026f27387fe0
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
