#!/bin/bash
cd ~/paw/moller
script=asym240
temp=${script}_temp.kumac
script=${script}.kumac
start=0
end=0

if [ "$#" -gt 0 ]; then
    start=$1
    end=$start
    echo "Start run $start"
fi
if [ "$#" -gt 1 ]; then
    end=$2
    echo "End run $end"
fi
for i in $(seq $start 1 $end);do
    echo "Run $i"
    cp $script $temp
    sed -i "s/run=0/run=$i/" $temp
    sed -i "s/\*shell .\/populate_runs_db.sh/ shell .\/populate_runs_db.sh $i/" $temp
    paw -b $temp
    #rm $temp
done
