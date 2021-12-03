FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:5268974eb6b041929d11f08b7bebc30d22af56ef7805703834c84aebfa437329
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
