FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:975b6e3a43e93ea254e0ea9eab49b1f4c6bd8d30c60def0416761bd36adf9acc
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
