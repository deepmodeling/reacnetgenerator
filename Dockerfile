FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:659a7ace1dc42fe929fe39a2ee729cefb2ccafe0b9fc8b04f0c3b575ab5d1ce8
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
