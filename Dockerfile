FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:48a57a493dd7fc63a60330d8bd96ffb2ef4093b0df1d413fb9b1ade4dac17f9d
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
