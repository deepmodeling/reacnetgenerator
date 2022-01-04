FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:062092f13b012cb6aeae3ed559afb866277298d3c140e039e4267dad15f88b40
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
