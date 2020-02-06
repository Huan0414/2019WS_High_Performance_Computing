#!/bin/sh
if [ `uname` = "Darwin" ]; then
BIN_DIR=
else
platform=$(uname -m)
if [ "$platform" = "x86_64" ]; then
BIN_DIR=bin64
export PKG_CONFIG_PATH='/home/ghaase/Lectures/Math2CPP/Codes/par/mgrid/'include/pkgconfig/lib64:$PKG_CONFIG_PATH
elif [ "$platform" = "amd64" ]; then
BIN_DIR=bin64
export PKG_CONFIG_PATH='/home/ghaase/Lectures/Math2CPP/Codes/par/mgrid/'include/pkgconfig/lib64:$PKG_CONFIG_PATH
else
BIN_DIR=bin32
export PKG_CONFIG_PATH='/home/ghaase/Lectures/Math2CPP/Codes/par/mgrid/'include/pkgconfig/lib32:$PKG_CONFIG_PATH
fi
fi
export PATH='/home/ghaase/Lectures/Math2CPP/Codes/par/mgrid/'$BIN_DIR:$PATH
if [ "$1" != "quiet" ]; then
echo "Copyright (C) 2009-2019 Intel Corporation. All rights reserved."
echo "Intel(R) Inspector 2019 (build 602103)"
fi
export INSPECTOR_2019_DIR='/home/ghaase/Lectures/Math2CPP/Codes/par/mgrid'
