FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:b4709f68cec496378ff17dcfa0458a14f3f955b314069a78166adcc0e89a5183
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
