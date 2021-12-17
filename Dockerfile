FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:89e630d9de7b439b90406e4793fcc0c423adb59ac7b02d7cccb5438b429684d6
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
