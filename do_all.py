#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 09:35:47 2018

@author: a123
"""

import os
import sys
import glob
import subprocess
import numpy as np

## First, extract the sac files and the resp files from the seed
## this can be done using "rdseed -R -z 1 -o 1 -d -f *.seed"
## and then "mkdir RESP_file;mv RESP* RESP_file";
## "cat RESP.*.*.*.* >> RESP.ALL" 
## then run this python file to process the data
## input dir where the sac files are to be stored
## e.g. "/Users/a123/Desktop/Repeating_Ear/Chile"
## and the dir where the sac files are stored temporarily
## e.g. "/Users/a123/Desktop/Repeating_Ear/Chile/unused"
os.putenv("SAC_DISPLAY_COPYRIGHT", '0');

if len(sys.argv) != 2:
    sys.exit("Usage: python %s dirname" %sys.argv[0]);

dir = sys.argv[1];
dir1 = sys.argv[2];

os.chdir(dir);
os.system('mkdir plot');
dir2 = dir1 + "/"
dir1 = dir1 + "/*.SAC";


p = subprocess.Popen(['sac'],stdin=subprocess.PIPE);
s = " ";
for fname in glob.glob(dir1):
    s += "r %s \n" %fname ;
    s += "p1 \n" ;
    s += "pause \n";
    s += "qdp off \n";
    s += "rmean \n" ;
    s += "rtrend \n";
    s += "trans from evalresp fname RESP_file/RESP.ALL to vel \n";
    s += "int \n";
    s += "p1 \n";
    s += "pause \n";
    s += "traveltime picks 0 phase P \n";
    s += "p1 \n"; 
    s += "pause \n";
    s += "p1 \n";
    s += "lh t0MARKER \n";
    s += "p1 \n"; 
    s += "pause \n"; 
    s += "setbb response (REPLY \"Hit 1 to do the polarity reversal \") \n";
    s += "if %response EQ 1 \n";
    s += "mul -1 \n";
    s += "p1 \n";
    s += "endif \n";
    s += "setbb response (REPLY \"Hit 1 to use this waveform for analysis \") \n";
    s += "if %response EQ 1 \n";
    s += " ppk BELL OFF \n";
    s += " setbb response1 (REPLY \"Hit 1 to write the picked times to header: ) \n";
    s += " if %response1 EQ 1 \n";
    s += "   setbb FILE1 (DELETE" + " " + dir2 + " " + fname + ")" + " " + "\n"
    s += "   write %FILE1 \n";
    s += " endif \n";
    s += "endif \n";
    s += "mv %s unused/already_checked \n" %fname; 

s += "q \n";
p.communicate(s.encode());
input("Press Enter to continue...")
## begin to the section of plotting 

os.chdir('plot');
os.system('cp ../*.SAC .');
## change the header 
reminder = "change the beginning time ..."
print(reminder)
p = subprocess.Popen(['sac'],stdin=subprocess.PIPE);
s = " ";
for fname in glob.glob("*.SAC"):
    s += "r %s \n" %fname;
    s += "ch B (&1,B& - &1,t0& + 1) \n";
    s += "ch t0 1 \n";
    s += "write over \n"; 
s += "q \n";
p.communicate(s.encode());

## do the resampling
delta = np.array([])
iii = 0; 
for dd in os.popen("saclst delta f *.SAC | awk '{print $2}' ").readlines():
    delta[iii] = float(dd.split()[0]);
    iii += 1;
sampling_rate = np.amax(delta, axis=0);

reminder = "now start to resample the data to %f/s" %sampling_rate;
print(reminder);

p = subprocess.Popen(['sac'],stdin=subprocess.PIPE);
s = " "; 
for fname in glob.glob("*.SAC"):
    s += "r %s \n" %fname;
    s += "interp delta %f \n" %sampling_rate;
    s += "write over \n";
s += "q \n";
p.communicate(s.encode());

### first pick the temp sacfile
input("Press Enter to continue...")
temp = input("pick the temp file... ")
t1mark = os.popen("saclst b f " + temp + " | awk '{print $2}' ").readline()
t1mark = t1mark.split()[0];
t0 = input("the onset of time window: ");
t0 = str(t0);
t1 = input("the end of time window: ");
t1 = str(t1);
## Before we run the xcorr_shift, please make sure it has been compiled, 
cmd = "xcorr_shift.o" + " " + "saclist" + " " + temp + " " + t0 + " " + t1 + " " + t1mark 
os.system(cmd)

os.system("saclst STLO STLA f *.SAC | awk '{print $2,$3}' >  sta.lst");
### run the gmt file...
### change the average time marker, lon, lat of the event
os.system('./plot.gmt');





    
