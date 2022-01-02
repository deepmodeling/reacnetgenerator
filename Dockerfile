FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:a83d15f55f15dac0ca34086b466575d3eced1222a21a34535401a049a9245c04
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
