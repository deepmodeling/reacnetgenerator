FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:a0fedf065ef4ad1c39f07df94ddd49c88634c5f27c65a3daea9b91e7e0de96e7
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
