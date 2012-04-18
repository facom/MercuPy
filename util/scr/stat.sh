#!/bin/bash
. sci2web/bin/sci2web.sm

#GET THE TOTAL SIMULATION TIME 
if [ -e "stdout.oxt" ];then
    time=$(cat stdout.oxt | grep Time: | tail -n 1 | awk -F':' '{print $2}' | awk '{print $1}')

    #CHECK PREVIOUS PHASES
    if [ -e "output/phase" ];then
	phase=$(cut -f 1 -d ' ' output/phase)
    else
	phase=0
    fi
    
    #CONTINUATION
    phase=$((phase%1))
    
    #COMPUTE THE FRACTION OF THE TOTAL TIME COMPLETED
    if [ 1 -gt 1 ];then
	elapsed=$(calc "100.0*1")
    else
	elapsed=100.0
    fi
    
    #ADVANCED STATUS
    stat=$(calc "($time+100.0*$phase)/$elapsed")
else
    stat=-1
fi

#REPORT STATUS
echo $stat
