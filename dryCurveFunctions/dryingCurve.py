# Generate consolidated dataframe from multiple CSVs
def concat(path, TrayW):
    all_files = glob.glob(path + "/*.csv")
    list = []

    #Reading CSV to dataframe
    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        df.drop(df.iloc[:, 12:19], inplace=True, axis=1)
        df.drop(df.iloc[:, 0:1], inplace=True, axis=1)
        df.drop(df.index[0], inplace=True)
        df.drop(df.index[0], inplace=True)
        list.append(df)

    frame = pd.concat(list, axis=0, ignore_index=True)
    frame.fillna(0, inplace=True)
    frame["W"] = frame["Sample_weight_g"] - TrayW
    frame["t"] = np.linspace(0,10*len(frame),len(frame)+1)[:-1] /3600
    
    return frame 


#frame["err"] = np.exp(-1*0.0214*frame["t"])*27.906
#frame["W"] = frame["Sample_weight_g"] - TrayW + frame["err"]
#print(min(frame["err"]), max(frame["err"]))

# Averaging for nset reading
def clean(par, nset):
    ret = []; k = 0; temp1 = 0; temp3=[]; temp2=[]; 
    
    while k < len(par):
        temp1 = 0; i = 0;
        while i <= nset:
            if (k >= len(par)): break
            temp2.append(par[k]); #print(par[k])
            i = i + 1; k = k + 1
        
        m = np.mean(temp2); #print('mean=', m)
        s = np.std(temp2);  #print('stdev=', s)
        
        for e in temp2:
                if (m - 3*s < e < m + 3*s) and (e > 0):
                    temp3.append(e)
        if len(temp3) > 0:
            temp1 = float("%.4f" % (sum(temp3)/len(temp3))); #print (temp1)
        else:
            temp1 = 0
        temp3=[]; temp2=[]
        ret.append(temp1); #print (ret)
        
        for i in range(len(ret)):
            if ret[i] == 0:
                ret[i] = m
    
    return ret
    del k, temp1, temp3, temp2, m, s, i, e, par
    
# Generate arrays for MC, MR
def genMCMR(wRed,exp, exDM):
    Water = np.subtract(wRed, exDM[exp]);         # print(min(Water),max(Water)) 
    MC = ( np.array(wRed) - exDM[exp] ) / wRed;   # print(min(MC),max(MC))
    MR = MC / MC[0];                              # print(min(MR),max(MR))
    return MC, MR
    
    # Plot data
def WeightVsTime(Time,Weight,xmin,xmax,xint,ymin,ymax,yint ):
    figure(num=None, figsize=(10, 8), dpi=300, facecolor='w', edgecolor='k')
    plt.scatter(Time, Weight,     linestyle='-', marker='o', color='k', linewidth=2)# , label='45 deg C | 0.2 m/s', s = 100) #, markersize='8'
    plt.xticks(np.arange(xmin,xmax,xint), fontsize=20); plt.xlabel('Time (hours)', fontsize=25)
    plt.yticks(np.arange(ymin,ymax,yint), fontsize=20); plt.ylabel('Weight (grams) ', fontsize=25)
    plt.grid(True); plt.ylim((ymin,ymax)); plt.xlim((xmin,xmax)); 
    #plt.legend(fontsize=20, loc='upper right', )
    plt.tight_layout(); plt.savefig('W_01.jpeg')
    # WeightVsTime(tRed,wRed,xmin,xmax,xint,ymin,ymax,yint)
    
# Plot data
def TempVsTime(Time,Temp,xmin,xmax,xint,ymin,ymax,yint):
    figure(num=None, figsize=(10, 8), dpi=300, facecolor='w', edgecolor='k')
    plt.scatter(Time, Temp,     linestyle='-', marker='o', color='r', linewidth=2)# , label='45 deg C | 0.2 m/s', s = 100) #, markersize='8'
    plt.xticks(np.arange(xmin,xmax,xint),  fontsize=20); plt.xlabel('Time (hours)', fontsize=25)
    plt.yticks(np.arange(ymin,ymax,yint), fontsize=20); plt.ylabel('Temperature (deg C) ', fontsize=25)
    plt.grid(True); plt.ylim((ymin,ymax)); plt.xlim((xmin,xmax));  
    #plt.legend(fontsize=20, loc='upper right', )
    plt.tight_layout(); plt.savefig('T_01.jpeg')

# Plot data
def RelHVsTime(Time,RelH,xmin,xmax,xint,ymin,ymax,yint):
    figure(num=None, figsize=(10, 8), dpi=300, facecolor='w', edgecolor='k')
    plt.scatter(Time, RelH,     linestyle='-', marker='o', color='g', linewidth=2)# , label='45 deg C | 0.2 m/s', s = 100) #, markersize='8'
    plt.xticks(np.arange(xmin,xmax,xint), fontsize=20); plt.xlabel('Time (hours)', fontsize=25)
    plt.yticks(np.arange(ymin,ymax,yint), fontsize=20); plt.ylabel('RH (deg C) ', fontsize=25)
    plt.grid(True); plt.ylim((ymin,ymax)); plt.xlim((xmin,xmax)); 
    #plt.legend(fontsize=20, loc='upper right', )
    plt.tight_layout(); plt.savefig('R_01.jpeg')
    
# Generate Moisture content Vs Time plot
def TimeVsMC(tRed, MC, MR, xmin,xmax,xint,ymin,ymax,yint):  
    figure(num=None, figsize=(10, 8), dpi=300, facecolor='w', edgecolor='k')   
    plt.scatter(tRed, MR, marker='+', color='k', label='MR')
    plt.scatter(tRed, MC, marker='^', color='r', label='MC')       
    plt.xticks(np.arange(xmin,xmax,xint), fontsize=20); plt.xlabel('Time (hours)', fontsize=25)
    plt.yticks(np.arange(ymin,ymax,yint), fontsize=20); plt.ylabel('(MC/MR) ', fontsize=25)
    plt.ylim((ymin,ymax)); plt.xlim((xmin,xmax)); plt.grid(True); 
    plt.legend(loc='upper right'); plt.tight_layout(); plt.savefig('MR_01.jpeg') #fontsize="small",        

import pandas as pd
import numpy as np
from numpy.polynomial import Polynomial
import glob
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pylab as pl
from sympy import S, symbols, printing
from matplotlib.lines import Line2D

# MarkerStyles, LineStyles, MarkerFills, Markercolors etc.
markerStyles = ['.','o','v','^','>','<','s','p','*','h','H','D','d','1','','+','X']
marker_style = dict(color='k', linestyle=' ', marker='s', markersize=15, markerfacecoloralt='b')
markerColors = ['b','g','r','c','m','y','k','w',]
T10colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
linestyle_str = [
     ('solid', 'solid'),      # Same as (0, ()) or '-'
     ('dotted', 'dotted'),    # Same as (0, (1, 1)) or '.'
     ('dashed', 'dashed'),    # Same as '--'
     ('dashdot', 'dashdot')]  # Same as '-.'
linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),

     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]