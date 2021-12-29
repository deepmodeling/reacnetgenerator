FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:1dae8d5044473cbde6fe1e80ee69aacd5f38d402b9a383eb7028653ffa5c7fc4
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
