FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:439127e006446ac8fe4986af2db3a4eed8bbe1a2caabc15017b4a6c1fee5413a
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
