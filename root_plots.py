import ROOT
from ROOT import TMath
import numpy as np


ttH = ROOT.TFile.Open("mc_341081.ttH125_gamgam.GamGam.root")
ggH = ROOT.TFile.Open("mc_343981.ggH125_gamgam.GamGam.root")
WWH = ROOT.TFile.Open("mc_345041.VBFH125_gamgam.GamGam.root")

data1 = ROOT.TFile.Open("data_A.GamGam.root")
data2 = ROOT.TFile.Open("data_B.GamGam.root")
data3 = ROOT.TFile.Open("data_C.GamGam.root")
data4 = ROOT.TFile.Open("data_D.GamGam.root")

Channels = {"ttH":ttH.Get("mini"),
            "ggH":ggH.Get("mini"),
            "WWH":WWH.Get("mini"), 
            'data1':data1.Get('mini'), 
            'data2':data2.Get('mini'),
            'data3':data3.Get('mini'),
            'data4':data4.Get('mini')
           }

for channel in Channels:
    print("Channel:",channel,"has",Channels[channel].GetEntries(),"entries")

OutputMap = {"ttH":[1,0,0,0],"ggH":[0,1,0,0],"WWH":[0,0,1,0], 
             'data1':[0,0,0,1], 'data2':[0,0,0,1], 'data3':[0,0,0,1], 'data4':[0,0,0,1]}

# Here we shall store the two photons & outputs 
dataset = []

# Prepare Data for NN
# Events are not filtered by LLT or HLT

for channel in Channels:
    print("Processing Channel ",channel)
    Channel = Channels[channel]
    counter = 0
    for event in Channel:
        #if (not event.trigP):
            #continue
        #if counter > 50000:
         #   break
        Photons = []
        if Channel.photon_n != 2:
            continue
        for j in range(Channel.photon_n):
            Momentum = ROOT.TLorentzVector()
            Momentum.SetPtEtaPhiE(Channel.photon_pt[j]/1000., Channel.photon_eta[j],Channel.photon_phi[j],Channel.photon_E[j]/1000.)
            Photons.append(Momentum)
        Photons.sort(key  = lambda p : -p.E())#sorts by energy most energetic goes first
        data = []
        for i in range(len(Photons)):
            #momentum of the photons go into nn
            data.append(Photons[i].E() )
            data.append(Photons[i].Px())
            data.append(Photons[i].Py())
            data.append(Photons[i].Pz())
        for vec in OutputMap[channel]:
            data.append(vec)
        dataset.append(data)
        counter += 1

canvas = ROOT.TCanvas('canvas', 'c', 800, 600) 

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
        
    #print(boundlow)
    #print(boundhigh)
    return hist

def photonpt(d):
	return np.sqrt(d[1]**2 + d[2]**2)

def eta(d):
    pt = np.sqrt(d[1]**2+d[2]**2)
    pz = d[3]
    return np.log(abs(np.tan((pt-pz)/(2*(pt+pz)))))

def invariantmass(d):
    return np.sqrt((d[0]+d[4])**2-(d[1]+d[5])**2-(d[2]+d[6])**2-(d[3]+d[7])**2)

def TInvariant(d):
    return np.sqrt(abs((d[0]-d[4])**2-(d[1]-d[5])**2-(d[2]-d[6])**2-(d[3]-d[7])**2))

def drawallhist(name, Lfunction):
    WWhist = drawhist('WWH_'+name, dataset, [lambda d: (True if d[10] == 1 else False)], lambda d: (Lfunction(d)))
    WWhist.Draw()
    ww = 'WWH'+name+'.jpeg'
    canvas.Print(str(ww))

    data1im = drawhist('dataset1'+name, dataset, [lambda d: (True if d[11] == 1 else False)], lambda d: (Lfunction(d)))
    data1im.Draw()
    ww = 'dataset1'+name+'.jpeg'
    canvas.Print(ww)

    data2im = drawhist('dataset2'+name, dataset, [lambda d: (True if d[11] == 1 else False)], lambda d: (Lfunction(d)))
    data2im.Draw()
    ww = 'dataset2'+name+'.jpeg'
    canvas.Print(ww)

    data3im = drawhist('dataset3'+name, dataset, [lambda d: (True if d[11] == 1 else False)], lambda d: (Lfunction(d)))
    data3im.Draw()
    ww = 'dataset3'+name+'.jpeg'
    canvas.Print(ww)

    data4im = drawhist('dataset4_'+name, dataset, [lambda d: (True if d[11] == 1 else False)], lambda d: (Lfunction(d)))
    data4im.Draw()
    ww = 'dataset4'+name+'.jpeg'
    canvas.Print(ww)

    ttHim = drawhist('tt'+name, dataset, [lambda d: (True if d[8] == 1 else False)], lambda d: (Lfunction(d)))
    ttHim.Draw()
    ww = 'ttH'+name+'.jpeg'
    canvas.Print(ww)

    ggHim = drawhist('gg'+name, dataset, [lambda d: (True if d[9] == 1 else False)], lambda d: (Lfunction(d)))
    ggHim.Draw()
    ww = 'ggH'+name+'.jpeg'
    canvas.Print(ww)

drawallhist('photon_pt', photonpt)
drawallhist('eta', eta)
drawallhist('invariant_mass', invariantmass)
drawallhist('Tinvariant', TInvariant)
