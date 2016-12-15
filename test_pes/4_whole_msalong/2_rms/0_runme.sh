#!/bin/bash
dir1 = ../1_fit


make
./rms.x *.fitting
rm rms.x

echo 'RMS vs Energy finished'
