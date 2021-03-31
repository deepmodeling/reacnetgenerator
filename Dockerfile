FROM python:3@sha256:d5d25f8ddcf983c0164bdcdc87b330d31417e2ce302dbd3e1d0e90fddf3ddff1
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
