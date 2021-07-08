FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:18a19b7c84e1a02ab58b167988695eaf915fc9ce7655346ebda2d3854dd47126
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
