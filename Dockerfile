FROM python:3
COPY . /reacnetgenerator
RUN pip install -–no-cache-dir /reacnetgenerator && \
    conda build purge-all && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
