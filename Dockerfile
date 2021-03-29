FROM python:3
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
