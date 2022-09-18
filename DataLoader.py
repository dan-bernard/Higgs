import ROOT
from ROOT import TMath
import random
import numpy as np

def LoadDataAndProcess(Channels, Outputs, Processing, MaxNumber = np.inf):
    dataset = []
    
    '''
        Prepare Data for NN
        Events are not filtered by LLT or HLT
    '''
    
    for channel in Channels:
        print("Processing Channel ",channel)
        Channel = Channels[channel]
        counter = 0
        for event in Channel:
            if counter > MaxNumber:
                break
            Photons = []
            if Channel.photon_n != 2:
                continue
            for j in range(Channel.photon_n):
                ''' 
                    Next we load the momentum vectors, the units within the root files 
                    are [MeV/c] for momentum and [MeV/c^2] for energy.
                    Here we use natural units where hbar = c = 1 and measure both 
                    energies and momenta in GeV, hence we need to divide by 1000 
                    the values coming from the root files.
                '''
                Momentum = ROOT.TLorentzVector()
                Momentum.SetPtEtaPhiE(Channel.photon_pt[j]/1000., Channel.photon_eta[j] ,Channel.photon_phi[j], Channel.photon_E[j]/1000.)
                Photons.append(Momentum)
            '''
                Sorting Function to always pass the photons to the NN 
                in order of decreasing energy.
            '''
            
            data = Processing(Photons)
            
            for vec in Outputs[channel]:
                '''
                    In addition to the input momenta and masses we
                    add the target outputs for the NN to the dataset
                '''
                data.append(vec)
            dataset.append(data)
            counter += 1
    
    '''
        We need to shuffle the data so that all the events coming 
        form a particular channel don't end up clumped together because 
        we later want to split the data into a training and testing set.
        We also save the data into a csv file for later quick access.
    '''
    
    random.shuffle(dataset)
    dataset = np.asarray(dataset)
    return dataset