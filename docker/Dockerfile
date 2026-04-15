FROM python:3.14

ENV PYTHONUNBUFFERED=1
EXPOSE 5000
WORKDIR /app

COPY ./requirements.txt .
RUN pip install -r requirements.txt

COPY . .

# Download and compile MAFFT
WORKDIR /app/mafft
RUN tar xzf mafft-7.525-with-extensions-src.tgz

# Modify the Makefile to fix bin path
RUN sed -i 's|^PREFIX = .*|PREFIX = /app/mafft/local|' mafft-7.525-with-extensions/core/Makefile && \
    sed -i 's|^BINDIR = .*|BINDIR = $(PREFIX)/bin|' mafft-7.525-with-extensions/core/Makefile && \
    sed -i 's|^LIBDIR = .*|LIBDIR = $(PREFIX)/lib|' mafft-7.525-with-extensions/core/Makefile

RUN cd mafft-7.525-with-extensions/core && \
    make clean && \
    make -j1 && \
    make install

# Ensure MAFFT_BINARIES is unset
ENV MAFFT_BINARIES=""

# Ensure the correct path is added to the environment
ENV PATH="/app/mafft/local/bin:$PATH"

WORKDIR /app
RUN chmod +x ./docker/entrypoint.sh
ENTRYPOINT ["./docker/entrypoint.sh"]