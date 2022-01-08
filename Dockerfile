FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:c25b61d4aa22abbc073e37e5fa25920a52c63a428587f2e72c473217d28577d2
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
