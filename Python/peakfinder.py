# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 13:00:34 2020

@author: Jérémie Marlhens
Title= peakfinder
"""
import numpy as np

from scipy.io import loadmat
x0 = loadmat('C:/Users/Grim.M/Ecole/ens/Master/S2/Stage/Code/cellular_cycle/x0.mat')
l=x0['x0'].tolist()
x=[item for sublist in l for item in sublist]
x0=x
"""
From Matlab file peakfinder.m: 
    
PEAKFINDER Noise tolerant fast peak finding algorithm
   INPUTS:
       x0 - A real vector from the maxima will be found (required)
       sel - The amount above surrounding data for a peak to be, 
           identified (default = (max(x0)-min(x0))/4). Larger values mean
           the algorithm is more selective in finding peaks.
       thresh - A threshold value which peaks must be larger than to be
           maxima or smaller than to be minima.
       extrema - 1 if maxima are desired, -1 if minima are desired
           (default = maxima, 1)
       includeEndpoints - If true the endpoints will be included as
           possible extrema otherwise they will not be included 
           (default = true)
       interpolate - If true quadratic interpolation will be performed 
           around each extrema to estimate the magnitude and the 
           position of the peak in terms of fractional indicies. Note that
           unlike the rest of this function interpolation assumes the 
           input is equally spaced. To recover the x_values of the input 
           rather than the fractional indicies you can do:
           peakX = x0 + (peakLoc - 1) * dx
           where x0 is the first x value and dx is the spacing of the
           vector. Output peakMag to recover interpolated magnitudes.
           See example 2 for more information.
           (default = false)

   OUTPUTS:
       peakLoc - The indicies of the identified peaks in x0
       peakMag - The magnitude of the identified peaks

   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
       are at least 1/4 the range of the data above surrounding data.

   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
       that are at least sel above surrounding data.

   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local 
       maxima that are at least sel above surrounding data and larger
       (smaller) than thresh if you are finding maxima (minima).

   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
       data if extrema > 0 and the minima of the data if extrema < 0

   [peakLoc] = peakfinder(x0,sel,thresh,extrema, includeEndpoints) 
       returns the endpoints as possible extrema if includeEndpoints is
       considered true in a boolean sense

   [peakLoc, peakMag] = peakfinder(x0,sel,thresh,extrema,interpolate) 
       returns the results of results of quadratic interpolate around each
       extrema if interpolate is considered to be true in a boolean sense

   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
       local maxima as well as the magnitudes of those maxima

   If called with no output the identified maxima will be plotted along
       with the input data.

   Note: If repeated values are found the first is identified as the peak

 Example 1:
 t = 0:.0001:10;
 x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
 x(1250:1255) = max(x);
 peakfinder(x)
 
 Example 2:
 ds = 100;   Downsample factor
 dt = .001;  Time step
 ds_dt = ds*dt;  Time delta after downsampling
 t0 = 1;
 t = t0:dt:5 + t0;
 x = 0.2-sin(0.01*2*pi*t)+3*cos(7/13*2*pi*t+.1)-2*cos((1+pi/10)*2*pi*t+0.2)-0.2*t;
 x(end) = min(x);
 x_ds = x(1:ds:end);  Downsample to test interpolation
 [minLoc, minMag] = peakfinder(x_ds, .8, 0, -1, false, true);
 minT = t0 + (minLoc - 1) * ds_dt;  Take into account 1 based indexing
 p = plot(t,x,'-',t(1:ds:end),x_ds,'o',minT,minMag,'rv'); 
 set(p(2:end), 'linewidth', 2);  Show the markers more clearly
 legend('Actual Data', 'Input Data', 'Estimated Peaks');
 Copyright Nathanael C. Yoder 2015 (nyoder@gmail.com)

"""
def peakfinder(x0, sel, thresh, extrema=1, includeEndpoints=0, interpolate=0): 
    
    if sel==[]:
        sel = (max(x0)-min(x0))/4;
    
        
    len0=len(x0)
    x0 = extrema*x0 #Make it so we are finding maxima regardless
    thresh = thresh*extrema  #Adjust threshold according to extrema.
    dx0 = np.diff(x0) #Find derivative
    ind = [idx+1 for idx, val in enumerate([a * b for (a, b) in zip(dx0[0:-2], dx0[1:-1])]) if val < 0] #dx0[dx0 == 0]= -eps; #This is so we find the first of repeated values
     #Find where the derivative changes sign
    
    if includeEndpoints :
        x = x0[0] + [x0[i] for i in ind] + x0[-1];
        ind = [0] + ind + [len0-1];
        minMag = min(x);
        leftMin = minMag;
    else :
        x = [x0[i] for i in ind];
        minMag = min(x);
        leftMin = min(x[1], x0[1]);
        
    
    length=len(x)
    
    if length > 2: #Function with peaks and valleys
    #Set initial parameters for loop
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
            ii = -1;
        else :
            ii = 0;
            
        #Preallocate max number of maxima
        maxPeaks = np.ceil(length/2);
        peakLoc = np.zeros((int(maxPeaks),1));
        peakMag = np.zeros((int(maxPeaks),1));
        cInd = 0;
        #Loop through extrema which should be peaks and then valleys
        while ii < length-1 :
            ii = ii+1; #This is a peak
            #Reset peak finding if we had a peak and the next peak is bigger than the last or the left min was small enough to reset.
            if foundPeak :
                tempMag = minMag;
                foundPeak = 0;

                
            #Found new peak that was lager than temp mag and selectivity larger than the minimum to its left.
            if x[ii] > tempMag and x[ii] > leftMin + sel :
                 tempLoc = ii;
                 tempMag = x[ii];
                
                
            #Make sure we don't iterate past the length of our vector
            if ii == length-1 :
                break
            #We assign the last point differently out of the loop end
                
            ii = ii+1; #Move onto the valley
        #Come down at least sel from peak
            if not foundPeak and tempMag > sel + x[ii] :
                foundPeak = 1; #We have found a peak
                leftMin = x[ii];
                peakLoc[cInd] = tempLoc; #Add peak to index
                peakMag[cInd] = tempMag;
                cInd = cInd+1;
            elif x[ii] < leftMin :#New left minima
                leftMin = x[ii];
                
                #Check end point
        if includeEndpoints :
            if x[-1] > tempMag and x[-1] > leftMin + sel :
                peakLoc[cInd] = length;
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
        [peakMags,xInd] = [max(x), np.argmax(x)];
        if includeEndpoints and peakMags > minMag + sel :
            peakInds = ind[xInd];
        else :
            peakMags = [];
            peakInds = [];

    
#Apply threshold value.  Since always finding maxima it will always be larger than the thresh.
    if 'thresh' in locals() :
        m= [peak>thresh for peak in peakMags];
        peakInds = [peakInds[i] for i, x in enumerate(m) if x];
        peakMags = [peakMags[i] for i, x in enumerate(m) if x];
        

    if interpolate and not not peakMags :
        middleMask = (np.array([x > 1 for x in peakInds]) & np.array([x < len0 for x in peakInds])+0);
        noEnds = [peakInds[i] for i, x in enumerate(middleMask) if x];
        
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


    #Plot if no output desired
    return peakMags, peakInds
    
    
    
    
    
    
    