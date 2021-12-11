FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:fc4b2e2901e3adf1efcbcfab030c11c6dde7e94b929ef66a97e7831eb3b6c3c9
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
