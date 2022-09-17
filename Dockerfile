FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:caed6dc95a09558d5f69baea35d5ccdcfce89c8667b7f10a969016622024a1b2
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
