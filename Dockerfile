FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:04f02dab1c526e33c59d83a3501b51b77235fcf844b8c2395db5282fbf2df4fb
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
