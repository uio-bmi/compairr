FROM alpine:latest
WORKDIR /opt/compairr
COPY Makefile .
COPY src ./src
COPY test ./test
RUN apk add --no-cache libstdc++ make g++ && \
    make clean && make && make test && make install && make clean && \
    apk del make g++
ENTRYPOINT ["/usr/local/bin/compairr"]
CMD ["--help"]
