FROM ubuntu:20.04 AS builder

RUN apt update && apt install -y --no-install-recommends\
  g++ \
  make \
  python3 \
  zlib1g-dev

WORKDIR /src

ADD http://code.enkre.net/bgen/tarball/release/v1.1.7 v1.1.7.tgz

RUN tar -xzf v1.1.7.tgz \
 && cd v1.1.7 \
 && python3 waf configure \
 && python3 waf

COPY . /src/regenie

WORKDIR /src/regenie

RUN make BGEN_PATH=/src/v1.1.7

FROM ubuntu:20.04

RUN apt update && apt install -y --no-install-recommends libgomp1 \
  && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/regenie/regenie /usr/local/bin

ENTRYPOINT ["/usr/local/bin/regenie"]
