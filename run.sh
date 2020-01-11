#!/bin/bash

for i in $(seq -f%04g 4 87)
do
    echo "./gsort ./rootfile/run$i.root ./gsortfile/gsort$i.root"
done
