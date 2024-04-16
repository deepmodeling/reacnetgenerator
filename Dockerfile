FROM python:3.12 AS compile-image
RUN python -m pip install uv
RUN python -m uv venv /opt/venv
# Make sure we use the virtualenv
ENV PATH="/opt/venv/bin:$PATH"
ENV VIRTUAL_ENV="/opt/venv"
# Install package
COPY . /reacnetgenerator
RUN uv pip install /reacnetgenerator && \
    reacnetgenerator -h

FROM python:3.12 AS build-image
COPY --from=compile-image /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
CMD ["/bin/bash" ]
