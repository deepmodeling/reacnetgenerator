FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:78669f9572203fa093c8bea3651a803bd223638549aa9576b0187d89b2bdd90a
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
