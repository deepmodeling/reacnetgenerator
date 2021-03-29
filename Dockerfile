FROM python:3
COPY . /reacnetgenerator
RUN apt-get update && \
    apt-get install -y --no-install-recommends yarn && \
	rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir /reacnetgenerator && \
	apt-get remove -y yarn && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
