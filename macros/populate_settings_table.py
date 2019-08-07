#!/usr/bin/env python3
import sqlite3, sys
import numpy as np
from datetime import datetime as dt
start_run = 0
quiet = True
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

#loop through the files
for r in range(start_run, end_run+1):
    try:
        file = open("/adaqfs/home/moller/daq/coda2/RunInfo/mollerrun_{}.set".format(r), 'r')
        if quiet == False:
            print("/adaqfs/home/moller/daq/coda2/RunInfo/mollerrun_{}.set".format(r))
    except:
        if quiet == False:
            print("Error opening /adaqfs/home/moller/daq/coda2/RunInfo/mollerrun_{}.set".format(r))
        continue



    entry = {}
    #Find run number
    run = 0
    for i in file:
        if "Run Number" in i:
            run = int(i.split()[3])
            if quiet == False:
                print("Run number: ",run)
            entry['run'] = int(run)
            break
        
    
    #Find run type
    run_type = ""
    for i in file:
        if "Type" in i:
            run_type = i[i.find(':')+1:-1].strip()
            if quiet == False:
                print("Run type: ",run_type)
            entry['run_type'] = str(run_type)
            break
    if run_type == 'moller_puls':
        print("Skipping pulser runs")
        continue
            
    #Find run start time
    run_start = 0
    for i in file:
        if "Date" in i:
            st = i[i.find(':')+1:-1].strip()
            run_start = dt.strptime(st,'%a %b %d %H:%M:%S %Z %Y')
            if quiet == False:
                print("Start time: ",run_start)
            entry['run_start'] = str(run_start)
            break
            
            
    #Find TD8000 discriminator thresholds
    found_group = False
    trig_thresh = []
    for i in file:
        if "Module: TD8000" in i or found_group:
            found_group = True
            if "THRESHOLD" in i:
                trig_thresh.append(float(i.split()[1]))
                trig_thresh.append(float(i.split()[2]))
                entry['trig_thresh_ch0'] = trig_thresh[0]
                entry['trig_thresh_ch1'] = trig_thresh[1]
                if quiet == False:
                    print("Discriminator thresholds: left=",trig_thresh[0]," mV,  right=",trig_thresh[1]," mV")
                break
                        
        
    #Find trigger configuration (Lecroy-2365)
    found_group = False
    trig_type = ''
    for i in file:
        if "Module: LeCroy-2365" in i or found_group:
            found_group = True
            if '7         4097     8192    1    not(not(0)+not(12)+(13))' in i:
                trig_type = 'left'
                entry['trig_type'] = str(trig_type)
                if quiet == False:
                    print('Trigger type: ', trig_type) 
                break
            if '7         4098     8192    1    not(not(1)+not(12)+(13))' in i:
                trig_type = 'right'
                entry['trig_type'] = str(trig_type)
                if quiet == False:
                    print('Trigger type: ', trig_type) 
                break
            if '7         4099     8192    1    not(not(0)+not(1)+not(12)+(13))' in i:
                trig_type = 'coinc'
                entry['trig_type'] = str(trig_type)
                if quiet == False:
                    print('Trigger type: ', trig_type) 
                break

    #If trigger type not found go back to start of file
    if trig_type == '':
        file.seek(0)
                            
                            
    #Load Epics information
    var_names ={'E_beam':'Beam energy, MeV Hall A',
                'E_inj':'Injector energy, MeV',
                'E_Slinac':'South linac energy, MeV',
                'E_Nlinac':'North linac energy, MeV',
                'n_pass': 'Passes Hall A',
                'bcm_avg': 'Beam Current Average',
                'unser': 'Current on Unser monitor',
                'bcm_us': 'Current on Upstream bcm',
                'bcm_ds': 'Current on Downstream bcm',
                'inj_bcm_tot': 'Injector Full Current Monitor 02',
                'inj_bcm_halla': 'Injector Current Monitor Hall A',
                'bpm01_X': 'Beam Position BPM01  X, mm',
                'bpm01_Y': 'Beam Position BPM01  Y, mm',
                'bpm04_X': 'Beam Position BPM04  X, mm',
                'bpm04_Y': 'Beam Position BPM04  Y, mm',
                'bpm04a_X': 'Beam Position BPM04A  X, mm',
                'bpm04a_Y': 'Beam Position BPM04A  Y, mm',
                'q1_cur': 'Quad Q1 (Amps)',
                'q2_cur': 'Quad Q2 (Amps)',
                'q3_cur': 'Quad Q3 (Amps)',
                'q4_cur': 'Quad Q4 (Amps)',
                'dip_cur': 'Dipole  (Amps)',
                'tgt_angle':'Target Rotary Position(V)',
                'tgt_angle_deg':'Rotary Position in deg',
                'tgt_lin_pos':'Target Linear Position(V)',
                'tgt_lin_pos_mm':'Linear Position in mm',
                'ihwp':'Laser 1/2 wave plate',
                'rhwp':'Rotating 1/2 wave plate',
                'vwien_angle':'VWien filter angle, deg',
                'sol_phi_fg':'Solenoids angle, deg',
                'hwien_angle':'HWien filter angle, deg',
                'hel_pattern':'Helicity Mode ON/OFF Random/Toggle',
                'hel_freq':'Helicity frequency',
                'hel_delay':'Helicity delay',
                't_settle':'MPS signal, usec',
                't_stable':'Helicity window, usec',
                'bpm02a_X': 'Beam Position BPM02A X, mm',
                'bpm02a_Y': 'Beam Position BPM02A Y, mm',
                'mol_mag_cur_set':'AM430 Current Setpoint (A)',
                'mol_mag_cur_meas':'AM430 Measured Current (A)',
                'mol_mag_v_meas':'AM430 Measured Voltage (A)',
                'mol_mag_field_meas':'AM430 Measured Field (T)',
                'mol_mag_ramp_state':'AMS430 Ramp State',
                'mol_cooler_temp':'Cryocooler Temperature (K)',
                'mol_mag_T2temp':'Magnet(T2) Temperature(K)',
                'mol_mag_lead1_temp':'Magnet Lead #1 (T6) Temperature',
                'mol_mag_lead2_temp':'Magnet Lead #2 (T7) Temperature',
                'det_hv_ch1':'Moller Detector measured HV ch 1',
                'det_hv_ch2':'Moller Detector measured HV ch 2',
                'det_hv_ch3':'Moller Detector measured HV ch 3',
                'det_hv_ch4':'Moller Detector measured HV ch 4',
                'det_hv_ch5':'Moller Detector measured HV ch 5',
                'det_hv_ch6':'Moller Detector measured HV ch 6',
                'det_hv_ch7':'Moller Detector measured HV ch 7',
                'det_hv_ch8':'Moller Detector measured HV ch 8',
                'det_ap_ch1':'Moller Detector measured HV Ap 1',
                'det_ap_ch2':'Moller Detector measured HV Ap 2',
                'det_ap_ch3':'Moller Detector measured HV Ap 3',
                'det_ap_ch4':'Moller Detector measured HV Ap 4',
                'det_ap_ch5':'Moller Detector measured HV Ap 5',
                'det_ap_ch6':'Moller Detector measured HV Ap 6',
                'det_ap_ch7':'Moller Detector measured HV Ap 7',
                'det_ap_ch8':'Moller Detector measured HV Ap 8'}
    epics_vals = {}
    for k, t in var_names.items():
        for i in file:
            if t in i:
                try:
                    #Fails if not a number
                    epics_vals[k] = float(i.split(':')[-1])
                    entry[k] = epics_vals[k]
                    if quiet == False:
                        print(list(epics_vals.keys())[-1]," ",list(epics_vals.values())[-1])
                except:
                    #Interpret as string
                    epics_vals[k] = str(i.split(':')[-1]).strip()
                    entry[k] = epics_vals[k]
                    if quiet == False:
                        print(list(epics_vals.keys())[-1]," ",list(epics_vals.values())[-1])
                break


                    
    #Find Moller coils power supply current
    mol_pow_sup_cur = 0
    for i in file:
        if "Hcoils current (Amps)" in i:
            mol_pow_sup_cur = float(i.split(':')[-1])
            entry['mol_pow_sup_cur'] = float(mol_pow_sup_cur)
            if quiet == False:
                print("Moller power supply current: ", mol_pow_sup_cur)
            break
            
            
    #Find run end time
    run_end = 0
    for i in file:
        if "End Run Time" in i:
            st = i[i.find(':')+2:-1]
            run_end = dt.strptime(st,'%a %b %d %H:%M:%S %Z %Y')
            entry['run_end'] = str(run_end)
            if quiet == False:
                print("End time: ", run_end)
            break
            

    #Determine length of run
    if run_end !=0 and run_start !=0 :
        entry['run_length'] =  (run_end - run_start).total_seconds()
        if quiet == False:
            print("Run duration: ", entry['run_length'])

    #Determine which target is being used
    step = 35.0 #distance between targets
    tgt1_pos = 61.0 #position in mm of first target
    tgt_pos =  entry['tgt_lin_pos_mm']
    tol = 3.0 #allowed variation
    if tgt_pos > tgt1_pos-tol and  tgt_pos < tgt1_pos+tol:
        entry['target'] = 1
    elif  tgt_pos > tgt1_pos+step-tol and  tgt_pos < tgt1_pos+step+tol:
        entry['target'] = 2
    elif  tgt_pos > tgt1_pos+2*step-tol and  tgt_pos < tgt1_pos+2*step+tol:
        entry['target'] = 3
    elif  tgt_pos > tgt1_pos+3*step-tol and  tgt_pos < tgt1_pos+3*step+tol:
        entry['target'] = 4
    else:
        entry['target'] = 0
    if quiet == False:
        print("Target number: ", entry['target'])


    #Make sql query
    stmt = "INSERT OR REPLACE INTO moller_settings("
    values = " VALUES("
    
    for k, t in entry.items():
        stmt = stmt + k +','
        values = values + '?,'
    stmt = stmt.rstrip(',')
    stmt = stmt + ')'
    values = values.rstrip(',')
    values = values +')'
    stmt = stmt + values
    if quiet == False:
        print(stmt)
        print(tuple(entry.values()))
    try:
        c.execute(stmt, tuple(entry.values()))
        conn.commit()
        if quiet == False:
            print("Successful insertion into moller_settings for run ", r)
    except:
        if quiet == False:
            print("FAILED insertion into moller_settings for run ", r)

conn.close()

