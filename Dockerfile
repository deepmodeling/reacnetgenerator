FROM nikolaik/python-nodejs:python3.9-nodejs14@sha256:299a588fbca99e907a3380d132487145a6032ba25983a6252179af90ed81b751
COPY . /reacnetgenerator
RUN pip install --no-cache-dir /reacnetgenerator && \
    reacnetgenerator -h
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["/bin/bash" ]
