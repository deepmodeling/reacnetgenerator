FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:1bdac451388c6442c0c92493d5341b59a8d3e2c4a354f59d1034e47b83b6def3
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
