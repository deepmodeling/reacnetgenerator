FROM python:3@sha256:cfa59e777e54fd12e9c90c9241ff19b5e9e07ce3db0f7e14763db5aceb0d3e28
COPY . /reacnetgenerator
RUN curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | apt-key add - && \
    echo "deb https://dl.yarnpkg.com/debian/ stable main" | tee /etc/apt/sources.list.d/yarn.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends nodejs yarn && \
	rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir /reacnetgenerator && \
	apt-get remove -y nodejs yarn && \
	apt-get autoremove -y && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
