FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:55a6e6a39f0317fc1e3c54151ba1153ccf3c5488ee0fb3b174334e6c10feab7c
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
