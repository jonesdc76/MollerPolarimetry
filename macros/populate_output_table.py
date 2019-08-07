#!/usr/bin/env python3
import sqlite3, sys
import numpy as np
start_run = 0
quiet =  True
if (len(sys.argv) > 1 ):
    start_run = int(sys.argv[1])
if  (len(sys.argv) > 2):
    end_run = int(sys.argv[2])
else: 
    end_run = start_run

#Connect to the existing database
try:
    conn=sqlite3.connect('/adaq1/data1/moller/MollerRunsDB.sql')
    c = conn.cursor()
except:
    print("Could not connect to database. Exiting")
    exit()

#names of variables in moller_output table
names = ('run','left_singles','right_singles','coinc','accid','bcm','clock','corrected_asym','corrected_asym_err','pol','pol_err','angle','analyzing_pow','target_pol','pol_left','pol_right','bcm_asym','bcm_asym_err')

#loop through the files
for r in range(start_run, end_run+1):
    try:
        file = open("/adaq1/data1/moller/.moller_{}.out".format(r),'r')
        if quiet == False:
            print("Opening /adaq1/data1/moller/.moller_{}.out".format(r))
    except:
        if quiet == False:
            print("Error opening file .moller_{}.out".format(r))
        continue

    for line in file:
        if (line.find(str(r)) == 0):
            entry = {names[i] : float(line.split()[i]) for i in range(len(names))}
            entry.pop('angle',None)
            if quiet == False:
                print(entry)
            try:
                c.execute("INSERT OR REPLACE INTO moller_output(run,left_singles,right_singles,coinc,accid,bcm,clock,corrected_asym,corrected_asym_err,pol,pol_err,analyzing_pow,target_pol,pol_left,pol_right,bcm_asym,bcm_asym_err) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", tuple(entry.values()))
                if quiet == False:
                    print("Successful insertion into moller_output for run ", r)
                break
            except:
                if quiet == False:
                    print("FAILED insertion into moller_output for run ", r)
                break

conn.commit()
conn.close()
