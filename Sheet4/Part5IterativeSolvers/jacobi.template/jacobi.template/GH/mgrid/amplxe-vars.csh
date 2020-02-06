#!/bin/csh
if ( `uname` == "Darwin" ) then
 set BIN_DIR=
else
  set platform = `uname -m`
  if ( ${platform} == "x86_64" ) then
 set BIN_DIR=bin64
  else if ( ${platform} == "amd64" ) then
 set BIN_DIR=bin64
  else
 set BIN_DIR=bin32
  endif
endif
set path = ("/home/ghaase/Lectures/Math2CPP/Codes/par/mgrid/${BIN_DIR}" $path)
if ( "$1" != "quiet" ) then
echo "Copyright (C) 2009-2019 Intel Corporation. All rights reserved."
echo "Intel(R) VTune(TM) Amplifier 2019 (build 602217)"
endif
setenv VTUNE_AMPLIFIER_2019_DIR "/home/ghaase/Lectures/Math2CPP/Codes/par/mgrid"
