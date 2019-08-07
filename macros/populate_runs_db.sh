#!/bin/bash
start=0
end=0
if [ "$#" -gt 0 ]; then
    start=$1
    end=$start
fi
if [ "$#" -gt 1 ]; then
    end=$2
fi

python3 populate_output_table.py $start $end
python3 populate_settings_table.py $start $end
