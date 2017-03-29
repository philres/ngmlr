FROM alpine

MAINTAINER Philipp Rescheneder <philipp.rescheneder@gmail.com>

ARG VERSION_ARG

# Install curl
RUN apk --update upgrade && \
    apk add git cmake curl ca-certificates && \
    update-ca-certificates && \
    rm -rf /var/cache/apk/*

# Get most recent ngmlr version from github
RUN echo $VERSION_ARG
#RUN curl -L `curl -s https://api.github.com/repos/philres/ngmlr/releases | grep browser_download_url | head -n 1 | cut -d '"' -f 4` | tar xz && cp ./ngmlr-*/ngmlr /bin/ && rm -rf ngmlr-*
RUN git clone https://github.com/philres/ngmlr.git && cd ngmlr && git checkout $VERSION_ARG && mkdir -p build && cd build && cmake .. && make && cp ../bin/ngmlr-*/ngmlr /bin/ && cd .. && rm -rf ngmlr 