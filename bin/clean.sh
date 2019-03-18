#!/bin/bash

export GOPATH=$(pwd) 
go clean src/main.go 
rm -f alloc
exit 0
