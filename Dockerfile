FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:cc31eaf8d387167737a11a4f8b7b5b3f425228b83e50ebd3286b6c2d76f43c02
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
