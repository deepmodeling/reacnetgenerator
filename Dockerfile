FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:875843f2f21e451f0919adc07efdff2ddc481c6620f02348efa16e1833dfd08f
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
