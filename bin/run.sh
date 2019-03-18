#!/bin/bash

echo "prepareing evn..."
export GOPATH=$(pwd)
echo "executing..."
go run ./src/main.go 
exit 0

