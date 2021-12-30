FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:aae4f6b02d86cf8d1d4350d2e3d186be8973b3fd939e1d2275a916d9436141a1
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
