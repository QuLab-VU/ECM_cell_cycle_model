# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 13:00:34 2020

@author: Jérémie Marlhens
Title= peakfinder
"""


import numpy as np
import sys
#
#x0 = cells['molecules'][7]['obsCycE']
#x0 = cells['molecules'][8]['obsCdc20']
#sel = 0.15*omega
#thresh = 0.15*omega
#extrema = 1
#includeEndpoints=1
#interpolate=0
#
#
#peakfinder(cells['molecules'][7]['obsCycE'], sel, thresh, extrema=1, includeEndpoints=1,interpolate=0)
#peakfinder(x0, sel, thresh, extrema=1, includeEndpoints=1,interpolate=0)
# 
# 
 
def peakfinder(x0, sel, thresh, extrema=1, includeEndpoints=0, interpolate=0): 
    
    len0 = len(x0)
        
    if sel==[]:
        sel = (max(x0)-min(x0))/4;
    
    
    x0 = extrema*x0 #Make it so we are finding maxima regardless
    thresh = thresh*extrema  #Adjust threshold according to extrema.
    dx0 = np.diff(x0) #Find derivative
    index=[idx for idx, val in enumerate(dx0) if val==0]
    for i in index :
        dx0[i]=-sys.float_info.epsilon
        
    ind = [idx+1 for idx, val in enumerate([a * b for (a, b) in zip(dx0[0:-2], dx0[1:-1])]) if val < 0] #dx0[dx0 == 0]= -eps; #This is so we find the first of repeated values
     #Find where the derivative changes sign
    
    if includeEndpoints :
        x = [x0[0]] + [x0[i] for i in ind] + [x0[-1]];
        ind = [0] + ind + [len0-1];
        minMag = min(x);
        leftMin = minMag;
        minLoc = x.index(minMag);
    else :
        x = [x0[i] for i in ind];
        if len(x)>0 :
            minMag = min(x);
            minLoc = x.index(minMag)
            x==minMag
            leftMin = min(x[0], x0[0]);
        
    length = len(x)
    
    if length > 2: #Function with peaks and valleys
    #Set initial parameters for loop
        tempLoc = minLoc
        tempMag = minMag;
        foundPeak = 0;
    
        if includeEndpoints :
            #Deal with first point a little differently since tacked it on
            #Calculate the sign of the derivative since we tacked the first 
            #point on it does not neccessarily alternate like the rest.
            signDx = np.sign(np.diff(x[0:3]));
            if signDx[0] <= 0 :#The first point is larger or equal to the second
                if signDx[0] == signDx[1] :#Want alternating signs
                    x[1] = [];
                    ind[1] = [];
                    length = length-1;
                    
            else : #First point is smaller than the second
                if signDx[0] == signDx[1] : #Want alternating signs
                    x[0] = [];
                    ind[0] = [];
                    length = length-1;
            
        
    
    
        #Skip the first point if it is smaller so we always start on a maxima
        if x[0] >= x[1] :
            kk = -1;
        else :
            kk = 0;
            
        #Preallocate max number of maxima
        maxPeaks = np.ceil(length/2);
        peakLoc = np.zeros((int(maxPeaks),1));
        peakMag = np.zeros((int(maxPeaks),1));
        cInd = 0;
        #Loop through extrema which should be peaks and then valleys
        while kk < length-1 :
            
            kk = kk+1; #This is a peak
            #Reset peak finding if we had a peak and the next peak is bigger than the last or the left min was small enough to reset.
            if foundPeak :
                tempMag = minMag;
                foundPeak = 0;
                
            #Found new peak that was lager than temp mag and selectivity larger than the minimum to its left.
            if x[kk] > tempMag and x[kk] > leftMin + sel :
                 tempLoc = kk;
                 tempMag = x[kk];
                
                
            #Make sure we don't iterate past the length of our vector
            if kk == length-1 :
                break
            #We assign the last point differently out of the loop end
                
            kk = kk+1; #Move onto the valley
        #Come down at least sel from peak
            if not foundPeak and tempMag > sel + x[kk] :
                foundPeak = 1; #We have found a peak
                leftMin = x[kk];
                peakLoc[cInd] = tempLoc; #Add peak to index
                peakMag[cInd] = tempMag;
                cInd = cInd+1;

            elif x[kk] < leftMin :#New left minima
                leftMin = x[kk];
                #Check end point
        
        
        if includeEndpoints :
            if x[-1] > tempMag and x[-1] > leftMin + sel :
                peakLoc[cInd] = length-1;
                peakMag[cInd] = x[-1];
                cInd = cInd + 1;
            elif not foundPeak and tempMag > minMag : #Check if we still need to add the last point
                peakLoc[cInd] = tempLoc;
                peakMag[cInd] = tempMag;
                cInd = cInd + 1;
        
        elif not foundPeak :
            if x[-1] > tempMag and x[-1] > leftMin + sel :
                peakLoc[cInd] = length;
                peakMag[cInd] = x[-1];
                cInd = cInd + 1;
            elif tempMag > min(x0[-1], x[-1]) + sel :
                peakLoc[cInd] = tempLoc;
                peakMag[cInd] = tempMag;
                cInd = cInd + 1;
                
        #Create output
        if cInd > 0 :
            peakInds = [ind[int(index)] for index in peakLoc[0:cInd]];
            peakMags = peakMag[0:cInd];
        else :
            peakInds = [];
            peakMags = [];
            
    else  :#This is a monotone function where an endpoint is the only peak
        if includeEndpoints :
            if length>0:
                [peakMags,xInd] = [[max(x)], np.argmax(x)];
            print(x0)
            if peakMags > minMag + sel :
                peakInds = [ind[xInd]];
            else :
                peakMags = [];
                peakInds = [];
        else :
            peakMags = [];
            peakInds = [];
    
#Apply threshold value.  Since always finding maxima it will always be larger than the thresh.
    if 'thresh' in locals() :
        m= [peak>thresh for peak in peakMags];
        peakInds = [peakInds[i] for i, xi in enumerate(m) if xi];
        peakMags = [peakMags[i] for i, xi in enumerate(m) if xi];
        

    if interpolate and not not peakMags :
        middleMask = (np.array([xi > 1 for xi in peakInds]) & np.array([xi < len0 for xi in peakInds])+0);
        noEnds = [peakInds[i] for i, x in enumerate(middleMask) if xi];
        
        magDiff = [(x0[i + 1] - x0[i - 1]) for i in noEnds];
        magSum = [(x0[i - 1] + x0[i + 1]  - 2 * x0[i]) for i in noEnds];
        magRatio = [x/y for x, y in zip(magDiff, magSum)];
        
        for i in middleMask.tolist() :
            peakInds[i] = peakInds[i] - magRatio/2;
            peakMags[i] = peakMags[i] - magRatio * magDiff/8;



    #Change sign of data if was finding minima
    if extrema < 0 :
        peakMags = -peakMags;
        x0 = -x0;


    return peakMags, peakInds
