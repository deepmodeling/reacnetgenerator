FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:db6611a3d16f1577353c766f93394467239a7f12cfd801c8825a25a5abb9e0aa
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
