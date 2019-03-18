#!/bin/bash

export GOPATH=$(pwd) 
rm -rf pkg
for item in $(ls $GOPATH/src); do [ -f "$GOPATH/src/$item" ] || rm -rf "$GOPATH/src/$item"; done
exit 0

