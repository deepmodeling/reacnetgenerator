FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:3b07a72cf71608d872f9b0f604b5ae5deaee710b562c96355b08735a8212a32c
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
