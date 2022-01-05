FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:4d8d6523bc2d74c9186456ab28dbd192d539fe9fbdd232089e3d98ba7e6033ae
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
