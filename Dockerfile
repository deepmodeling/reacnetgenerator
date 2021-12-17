FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:065863d4805f5c1311f50e8e19ae285015a82bfed74690ba4c77df3e53fea6ac
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
