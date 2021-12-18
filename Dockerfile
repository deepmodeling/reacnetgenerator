FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:c51319f60c82014123b5131ee5f609aa2d2ce22a2aa1e16475ab0dbead72ba3d
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
