FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:ada515fcd864804711aa5a183bb1781e4c760ca13204d1b52d71f3ae81c7f854
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
