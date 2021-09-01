FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:2b480f2a0a98fb781b27ade4b3ef22c85741a8e0fba78506747584cb0325c098
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
