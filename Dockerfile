# Thanks to Nathan Weeks for sharing this Dockerfile
# minor edits were made
# Filename: Dockerfile

# make this global 
ARG LIB_INSTALL
ARG LIB_INSTALL2


FROM ubuntu:18.04 AS builder

ARG BOOST_IO
ARG LIB_INSTALL
ARG STATIC

WORKDIR /src

ADD http://code.enkre.net/bgen/tarball/release/v1.1.7 v1.1.7.tgz

RUN apt-get update && apt-get install -y --no-install-recommends \
      g++ \
      make \
      python3 \
      zlib1g-dev \
      $LIB_INSTALL \
      && tar -xzf v1.1.7.tgz \
      && rm v1.1.7.tgz \
      && cd v1.1.7 \
      && python3 waf configure \
      && python3 waf

COPY . /src/regenie

WORKDIR /src/regenie

RUN make BGEN_PATH=/src/v1.1.7 HAS_BOOST_IOSTREAM=$BOOST_IO STATIC=$STATIC

FROM ubuntu:18.04
ARG LIB_INSTALL2

RUN apt-get update && apt-get install -y --no-install-recommends \
      libgomp1 $LIB_INSTALL2 \
      && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/regenie/regenie /usr/local/bin

# Avoid this to keep image more for general usage
# ENTRYPOINT ["/usr/local/bin/regenie"]

