FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:cb95117b1df47dc418889c4b82d1f970610cbe4a4b05f4c83580fb0960ca9ade
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
