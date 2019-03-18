#!/bin/bash

export GOPATH=$(pwd)
go build -o alloc -v -work -x src/main.go
exit 0

