FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:3083ab5d9d93bad3cd609470e7feb5c6ccc3df896b7e1aad818fdc458273715e
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
