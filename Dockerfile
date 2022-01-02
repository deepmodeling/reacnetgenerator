FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:9449cabcca8037b7a979e8b91e9e72eb988ef2306400a79a1c67b9ccf7be9524
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
