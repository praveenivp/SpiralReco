#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Flikering radial checkerboard pattern
Trigger should send 'q' or 'Q' or 'w'
press x to quit 
log file can be found in %TEMP%\\fmri_radialChecker_YYYYMMDD-HHMMSS.txt
"""
from __future__ import division

from psychopy import visual, event, core
import numpy as np
import logging,time,os
from random import randint

def  deg2pix2(degrees,ScreenResolution):
    # deg2pix2 converts visual angle in degress to fration of screen dimension 
    # pixels2d*ScreenResolution gives the actual number of pixels to cover the visual angle. 
    viewdist                                      = 1320.   # viewing distance (mm)
    screenh                                       = 262.00 # q height (mm)
    pixels = np.tan((degrees*np.pi/180))*viewdist/screenh
    aspect_ratio=ScreenResolution[0]/ScreenResolution[1]
    pixels2d= (pixels/aspect_ratio,pixels)
    return pixels2d

# All Parameters
screenResolution=(1920.,1200.)
# sequence timing: Now the Stimuli will be blank for #Pause_seg triggers and checkerbwwoard board will be displayed for #Stimuli_seg
Rest_trig=10 #triggers
Stimuli_trig=10 #triggers
total_trig=Stimuli_trig+Rest_trig
#stimuli parameters
contrast=0.8
StimSize=deg2pix2(20,screenResolution) # 20 deg to screen fraction
mask=np.array([1,1,1,1]) #4/4 of the circle is full
rotationRate = 0.00  # revs per sec
flashPeriod = 2/10 # seconds for one B-W cycle (ie 10 Hz)
Debug=False # print some helpful texts during debugging

#open a Log file and log the stimulus settings 
Fn = "fmri_radialChecker_"+ time.strftime("%Y%m%d-%H%M%S") +".txt"
tempdir=os.getenv("TEMP")
Fn=os.path.join(tempdir,Fn)
print("Logging to file : %s" % Fn)
LogFile=open(Fn, 'w')
LogFile.write("#Radial Checkerboard\n")
LogFile.write(f"Stimulus Paramters : \nRest/Stimulation Volumes = {Rest_trig:d}/{Stimuli_trig:d} Triggers \nFlash frequency= {2/flashPeriod:.4f} Hz \nContrast={contrast:.2f}\n" )
LogFile.write(('Mask='+''.join(['%d,']*mask.size))%tuple(mask)+'\n')
LogFile.write(('StimSize='+''.join(['%.4f,']*len(StimSize))+'\n\n\n')%(StimSize))
LogFile.write(f"{'Triggers':<10s}|{'Time(s)':<10s} |{'Stimulus':<10s}  \n")
LogFile.write('------------------------------\n')
LogFile.flush()

# window setup
win = visual.Window(screenResolution, allowGUI=False,screen=1, fullscr=True,color=-0.9)
globalClock = core.Clock()

#some text elements for debugging
trigtext = visual.TextStim(win, pos=(0.0, 0.9), text='Number of Trig: 0')
message = visual.TextStim(win, pos=(0.0, -0.8), text='Ready for Trigger(''q'' or ''w''). \n Hit X to quit' )
trigtext.autoDraw=Debug
message.autoDraw=Debug

# Make two checkerboard (in opposite contrast) and alternate them for flashing
checkerboard = visual.RadialStim(win, tex='sqrXsqr', color=contrast, size=StimSize, 
    visibleWedge=[0, 360], radialCycles=4, angularCycles=8, interpolate=False,
    autoLog=False, mask=mask)  # this stim changes too much for autologging to be useful
checkerboard_inv = visual.RadialStim(win, tex='sqrXsqr', color=-1*contrast, size=StimSize,
    visibleWedge=[0, 360], radialCycles=4, angularCycles=8, interpolate=False,
    autoLog=False,mask=mask)  # this stim changes too much for autologging to be useful
# fixation mark and blanck screen
fixation_point = visual.TextStim(win, pos=(0.0, 0), text='+', color=(0.43, 0.7, 0.28), colorSpace='rgb')
fixation_point.autoDraw=True
fixation = visual.GratingStim(win=win, size=0, pos=[0,0], sf=0, rgb=0)
# Define a list of the color values
Colors = [(255, 0, 0), (0, 255, 0), (0, 0, 255), (255, 255, 0), (255, 0, 255),(0, 255, 255), (255, 128, 0), (255, 192, 203), (139, 69, 19), (128, 128, 128)];
ColorChangePeroid=5 #seconds
ColorCount=0;


# run()
t = 0
nTrig=-1
State=0 # rest 1: stimulus 
ColCount=0;
print("ready for trigger\n")
print(f"{'Triggers':<10s}|{'Time(s)':<10s} |{'Stimulus':<10s}  \n")
timer = core.CountdownTimer(ColorChangePeroid)
while not event.getKeys('x'):
    t = globalClock.getTime()
    fixation_point.draw()
    if (event.getKeys(keyList=['q','w','Q']) ):
        nTrig=nTrig+1
        State=((nTrig % total_trig) >= Rest_trig )
        trigtext.text='Number of Trig: %d' % (nTrig+1)
        print(f"{nTrig+1:<10d} {t:<10.4f} {State:<10d}")  
        LogFile.write(f"{nTrig+1:<10d} {t:<10.4f} {State:<10d} \n")#('%d %.4f  %d \n' % (nTrig+1,t,State))
        LogFile.flush()
    if(timer.getTime()<0 and nTrig>-1):
        CurrCol=randint(0, 9);
        fixation_point.setColor(Colors[CurrCol],'rgb255')
        if(CurrCol==0): # let's count number of red crosses
            ColCount=ColCount+1
        timer.reset()
    if(State==1  and nTrig>0):
        if t % flashPeriod < (flashPeriod/2.0):  # more accurate to count frames
            stim = checkerboard
        else:
            stim = checkerboard_inv
        stim.ori = t * rotationRate * 360.0  # set new rotation
    else:
        stim=fixation
    stim.draw()
    win.flip()

print(f"\nNumber of red plus is {ColCount:<10d}")  
LogFile.close()
win.close()
core.quit()
