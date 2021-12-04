FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:2887f257ce9a70ff408fa798e0619c1d3d23dc7a03733a62d82e07ff8ad57020
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
