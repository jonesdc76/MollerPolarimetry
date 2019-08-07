#!/bin/bash

sqlite3 MollerRunsDB.sql <<EOF

CREATE TABLE IF NOT EXISTS moller_output(
run INTEGER PRIMARY KEY,           -- run number
left_singles INTEGER DEFAULT 0,    -- average left detector singles rate (per helicity)
right_singles INTEGER DEFAULT 0,   -- average right detector singles rate (per helicity)
coinc INTEGER DEFAULT 0,           -- average coincidence rate (per helicity)
accid INTEGER DEFAULT 0,           -- average accidental rate (per helicity)
bcm INTEGER DEFAULT 0,             -- average bcm scalar counts per helicity 
clock INTEGER DEFAULT 0,           -- average clock counts per helicity
corrected_asym REAL DEFAULT 0,     -- corrected asymmetry
corrected_asym_err REAL DEFAULT 0, -- corrected asymmetry error
pol REAL DEFAULT 0,                -- calculated polarization
pol_err REAL DEFAULT 0,            -- calculated polarization error
analyzing_pow REAL DEFAULT 0,      -- effective analyzing power Azz
target_pol REAL DEFAULT 0,         -- target polarization
pol_left REAL DEFAULT 0,           -- left detector only asymmetry
pol_right REAL DEFAULT 0,          -- right detector only asymmetry
bcm_asym REAL DEFAULT 0,           -- bcm asymmetry
bcm_asym_err REAL DEFAULT 0        -- bcm asymmetry error
);

CREATE TABLE IF NOT EXISTS moller_settings(
run INTEGER PRIMARY KEY,           -- run number
run_type TEXT DEFAULT '',          -- run type eg. beam_pol
run_start TEXT DEFAULT '',         -- run start date and time
run_end TEXT DEFAULT '',           -- run end date and time
run_length REAL DEFAULT 0,         -- run duration (seconds)
trig_thresh_ch0 REAL DEFAULT 0,    -- ch 0 (left) detector discriminator threshold (mV)
trig_thresh_ch1 REAL DEFAULT 0,    -- ch 1 (right) detector discriminator threshold (mV)
ihwp_in INTEGER DEFAULT -1,        -- insertable half-wave plate IN=1; OUT=0
target INTEGER DEFAULT 0,          -- target 0=not on any target; 1=Cu 11um; 2=Fe 10um; 3=Fe 4um; 4=Fe 1um;
E_beam REAL DEFAULT 0,             -- Beam energy, MeV Hall A
E_inj REAL DEFAULT 0,              -- Injector energy, MeV
E_Slinac REAL DEFAULT 0,           -- South linac energy, MeV
E_Nlinac REAL DEFAULT 0,           -- North linac energy, MeV
n_pass REAL DEFAULT 0,             -- Passes Hall A
bcm_avg REAL DEFAULT 0,            -- Beam Current Average
unser REAL DEFAULT 0,              -- Current on Unser monitor
bcm_us REAL DEFAULT 0,             -- Current on Upstream bcm
bcm_ds REAL DEFAULT 0,             -- Current on Downstream bcm
inj_bcm_tot REAL DEFAULT 0,        -- Injector Full Current Monitor 02
inj_bcm_halla REAL DEFAULT 0,      -- Injector Current Monitor Hall A
bpm01_X REAL DEFAULT 0,            -- Beam Position BPM01  X, mm
bpm01_Y REAL DEFAULT 0,            -- Beam Position BPM01  Y, mm
bpm04_X REAL DEFAULT 0,            -- Beam Position BPM04  X, mm
bpm04_Y REAL DEFAULT 0,            -- Beam Position BPM04  Y, mm
bpm04a_X REAL DEFAULT 0,           -- Beam Position BPM04A  X, mm
bpm04a_Y REAL DEFAULT 0,           -- Beam Position BPM04A  Y, mm
q1_cur REAL DEFAULT 0,             -- Quad Q1 (Amps)
q2_cur REAL DEFAULT 0,             -- Quad Q2 (Amps)
q3_cur REAL DEFAULT 0,             -- Quad Q3 (Amps)
q4_cur REAL DEFAULT 0,             -- Quad Q4 (Amps)
dip_cur REAL DEFAULT 0,            -- Dipole  (Amps)
tgt_angle REAL DEFAULT 0,          -- Target Rotary Position(V)
tgt_angle_deg REAL DEFAULT 0,      -- Rotary Position in deg
tgt_lin_pos REAL DEFAULT 0,        -- Target Linear Position(V)
tgt_lin_pos_mm REAL DEFAULT 0,     -- Linear Position in mm
ihwp TEXT DEFAULT '',              -- Laser 1/2 wave plate
rhwp REAL DEFAULT 0,               -- Rotating 1/2 wave plate
vwien_angle REAL DEFAULT 0,        -- VWien filter angle, deg
sol_phi_fg REAL DEFAULT 0,         -- Solenoids angle, deg
hwien_angle REAL DEFAULT 0,        -- HWien filter angle, deg
hel_pattern TEXT DEFAULT '',       -- Helicity mode/pattern
hel_freq REAL DEFAULT 0,           -- Helicity frequency
hel_delay TEXT DEFAULT '',         -- Helicity delay
t_settle REAL DEFAULT 0,           -- MPS signal, usec
t_stable REAL DEFAULT 0,           -- Helicity window, usec
bpm02a_X REAL DEFAULT 0,           -- Beam Position BPM02A X, mm
bpm02a_Y REAL DEFAULT 0,           -- Beam Position BPM02A Y, mm
mol_mag_cur_set REAL DEFAULT 0,    -- AM430 Current Setpoint (A)
mol_mag_cur_meas REAL DEFAULT 0,   -- AM430 Measured Current (A)
mol_mag_v_meas REAL DEFAULT 0,     -- AM430 Measured Voltage (A)
mol_mag_field_meas REAL DEFAULT 0, -- AM430 Measured Field (T)
mol_mag_ramp_state REAL DEFAULT 0, -- AMS430 Ramp State
mol_cooler_temp REAL DEFAULT 0,    -- Cryocooler Temperature (K)
mol_mag_T2temp REAL DEFAULT 0,     -- Magnet(T2) Temperature(K)
mol_mag_lead1_temp REAL DEFAULT 0, -- Magnet Lead #1 (T6) Temperature
mol_mag_lead2_temp REAL DEFAULT 0, -- Magnet Lead #2 (T7) Temperature
det_hv_ch1 REAL DEFAULT 0,         -- Moller Detector measured HV ch 1
det_hv_ch2 REAL DEFAULT 0,         -- Moller Detector measured HV ch 2
det_hv_ch3 REAL DEFAULT 0,         -- Moller Detector measured HV ch 3
det_hv_ch4 REAL DEFAULT 0,         -- Moller Detector measured HV ch 4
det_hv_ch5 REAL DEFAULT 0,         -- Moller Detector measured HV ch 5
det_hv_ch6 REAL DEFAULT 0,         -- Moller Detector measured HV ch 6
det_hv_ch7 REAL DEFAULT 0,         -- Moller Detector measured HV ch 7
det_hv_ch8 REAL DEFAULT 0,         -- Moller Detector measured HV ch 8
det_ap_ch1 REAL DEFAULT 0,         -- Moller Detector measured HV Ap 1
det_ap_ch2 REAL DEFAULT 0,         -- Moller Detector measured HV Ap 2
det_ap_ch3 REAL DEFAULT 0,         -- Moller Detector measured HV Ap 3
det_ap_ch4 REAL DEFAULT 0,         -- Moller Detector measured HV Ap 4
det_ap_ch5 REAL DEFAULT 0,         -- Moller Detector measured HV Ap 5
det_ap_ch6 REAL DEFAULT 0,         -- Moller Detector measured HV Ap 6
det_ap_ch7 REAL DEFAULT 0,         -- Moller Detector measured HV Ap 7
det_ap_ch8 REAL DEFAULT 0,         -- Moller Detector measured HV Ap 8'
trig_thresh REAL DEFAULT 0,        -- trigger threshold in mV
trig_type TEXT DEFAULT '',         -- trigger type: 'left', 'right' or 'coinc'
mol_pow_sup_cur REAL DEFAULT 0,    -- Moller power supply current (A)
FOREIGN KEY(run) REFERENCES moller_output(run)
);

CREATE TABLE IF NOT EXISTS moller_quality(
run INTEGER PRIMARY KEY,           -- run number
task TEXT DEFAULT '',              -- 'pol'=polarization, 'qsc'=quad scan, 'dsc'=dipole scan, 'thr'=threshold study, 'trg'=trigger study, 'hv'=high voltage study, 'ped'=pedestal, 'oth'=other
quality TEXT DEFAULT '',           -- 'good', 'suspect', 'bad'
comment TEXT DEFAULT '',           -- user input comment about run
FOREIGN KEY(run) REFERENCES moller_output(run)
);

EOF
