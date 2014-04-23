#!/bin/sh
# sh cmd.sh

for i in `ls -l /home/incompact3dNew/channel|grep f90|awk '{print $NF}'`;
do 
    file=/home/incompact3dNew/channel/$i;
#    echo $file;
     ./split $file
done;
