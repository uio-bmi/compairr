FROM alpine:3.13
WORKDIR /opt/compairr
COPY . .
RUN apk add --no-cache libstdc++ make g++ && \
    make clean && make && make test && make install && make clean && \
    apk del make g++
ENTRYPOINT ["/usr/local/bin/compairr"]
CMD ["--help"]
