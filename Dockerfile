FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:124f649a73801d4e1bd02088f218e747b15d55c5a7c4a828910df09bae2ed341
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
