FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:7d49857db799a58bc561044f40ab575e56854604157e77308a21ce7a653aeaba
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
