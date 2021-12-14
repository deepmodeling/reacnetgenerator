FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:21e3c32b964f9502b62aba767af6372da73371455fcd905b49012898eab0ca09
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
