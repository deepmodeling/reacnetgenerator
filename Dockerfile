FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:fb08cb56d03af6080c7e8a08dada0737a9cfa01c6d56063752a40110d533c1fd
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
