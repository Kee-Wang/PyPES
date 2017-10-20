#!/bin/bash
qstat -n >> /home/kee/github/PyPES/workspace/que/nodetemp
python2 /home/kee/github/PyPES/workspace/que/que_tool.py
rm /home/kee/github/PyPES/workspace/que/nodetemp
