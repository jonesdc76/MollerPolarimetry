#!/usr/bin/env python3
import sqlite3
import numpy as np
start_run = sys.argv[1]
end_run = sys.argv[2]
try:
    conn=sqlite3.connect('MollerRunsDB.sql')
    c = conn.cursor()
except:
    print("Could not connect to database. Exiting")
    return

#names of variables in moller_output table
names = ('run','left_singles','right_singles','coinc','accid','bcm','clock','corrected_asym','corrected_asym_err','pol','pol_err','analyzing_pow','target_pol','pol_left','pol_right','bcm_asym','bcm_asym_err')

#loop through the files
for r in range(start_run, end_run+1):
    try:
        file = open("/adaq1/data1/moller/.moller_{}.out".format(r))
    except:
        continue

    for line in file:
        if(int(line.split()[0]) == r):
            entry = {names[i] : float(i) for i in line.split()}
            try:
                c.execute("INSERT OR REPLACE INTO moller_output(run,left_singles,right_singles,coinc,accid,bcm,clock,corrected_asym,corrected_asym_err,pol,pol_err,analyzing_pow,target_pol,pol_left,pol_right,bcm_asym,bcm_asym_err) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", tuple(entry.values()))
            except:
                print("Insertion into moller_output failed for run ", r)

conn.commit()
conn.close()

 #   run, left, right, coin, acc, bcm, clock, corAsym, corAsymErr, pol, polErr, angl,  anPow, tgtPol, polL, polR, Abcm, AbcmErr, coil, factor = np.loadtxt(f'runlist.dat',unpack=True)
