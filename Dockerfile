FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:9b51cb43e457991e6ab83a3ec2cc00f2152a06b95def7a019567585250b0144a
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
