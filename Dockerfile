FROM python:3@sha256:f611a8d88dfec4e50f72fd3d80d082c778f56896243619e8d4eed3d719891ffc
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
