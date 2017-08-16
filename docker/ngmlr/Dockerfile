FROM alpine:3.3

MAINTAINER Philipp Rescheneder <philipp.rescheneder@gmail.com>

ARG VERSION_ARG

# Get most recent ngmlr version from github
RUN apk --update upgrade && \
    apk add build-base gcc abuild binutils binutils-doc gcc-doc zlib-dev git cmake curl ca-certificates && \
    update-ca-certificates && \
    git clone https://github.com/philres/ngmlr.git && cd ngmlr && git checkout $VERSION_ARG && mkdir -p build && cd build && cmake .. && make && cp ../bin/ngmlr-*/ngmlr /bin/ && cd .. && rm -rf ngmlr && \
    apk del build-base abuild binutils binutils-doc gcc-doc zlib-dev git cmake curl ca-certificates && \
    rm -rf /var/cache/apk/*
