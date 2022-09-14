FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:bc66bc69d782894ac3cb2e5bf0ff75348ebe27581760ceca3c3c1d5bac63ba31
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
