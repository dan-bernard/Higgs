import ROOT
import numpy as np

def prob(arr):
    return (arr[0] + arr[1] + arr[2])/sum(arr)

def drawhist(name, Data, filters, function):
    boundlow = 0
    boundhigh = 0
    
    for d in Data:
        skip = False
        for f in filters:
            if f(d) == False:
                skip = True
        if skip:
            continue
        g = function(d)
        if boundlow > g:
            boundlow = g
        if boundhigh < g:
            boundhigh = g
            
    
    hist = ROOT.TH1F(name, name, int(np.sqrt(len(Data))), boundlow, boundhigh)
    
    for d in Data:
        skip = False
        for f in filters:
            if f(d) == False:
                skip = True
        if skip:
            continue
        g = function(d)
        hist.Fill(g)
        
    print(boundlow)
    print(boundhigh)
    return hist