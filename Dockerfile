FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:2bdb590331b5be33d20430fcf6feeb808790d0bda50576ecae719adfa03d2637
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
