FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:993408995ded1b422d933046a774ad752a7b5880ce29368d64c49e0792d1528e
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
