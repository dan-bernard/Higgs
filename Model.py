# Turn off GPU Support 
import os
using_gpu = False
if(not using_gpu):
    os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

import tensorflow as tf
from tensorflow import keras
from keras import layers

import matplotlib.pyplot as plt
import random


class Model:

    def __init__(self,name,path,Dataset):
        ilayer = layers.Dense(units = 7, input_dim = len(Dataset[0,:-4]), activation = "relu")
        layer1 = layers.Dense(units = 8, activation = "relu")
        layer2 = layers.Dense(units = 9, activation = "relu")
        olayer = layers.Dense(units = 4, activation = "sigmoid")

        model = keras.Sequential([ilayer,layer1,layer2,olayer])
    
        model.compile(
        optimizer = keras.optimizers.SGD(0.001),
        loss="binary_crossentropy"
        )

        self.data  = Dataset
        self.name  = name
        self.model = model
        self.path  = path
        self.saver = tf.keras.callbacks.ModelCheckpoint(filepath = self.path, save_weights_only = True)

    def train(self,fraction,epochs,verbose=False):
        '''
        
            Uses the first fraction of the entire dataset as training
            @params:
                fraction: (double) fraction of the dataset to be trained on
        
        '''

        tf.keras.backend.clear_session()

        point = int(len(self.data)*fraction)

        self.xtrain = self.data[:point,:-4]
        self.ytrain = self.data[:point,-4:]

        self.xtest  = self.data[point:,:-4]
        self.ytest  = self.data[point:,-4:]

        history = self.model.fit(x = self.xtrain, y=self.ytrain, epochs = epochs, verbose = True, callbacks = [self.saver])
        plt.plot(history.history["loss"])
        plt.show()

    
    def setSavePath(self,path):
        self.path  = path
        self.saver = tf.keras.callbacks.ModelCheckpoint(filepath =self.path, save_weights_only = True)

    
    def loadFromPath(self,path):
        self.model.load_weights(path)

    def analyze(self):

        def process(prediction): #converts from % prediction to [0,0,1,0]
            for i in range(len(prediction)):
                Maxj = 0
                for j in range(len(prediction[i])):
                    if prediction[i][j] > prediction[i][Maxj]:
                        Maxj = j       
                prediction[i] = [0 if prediction[i][Maxj] > prediction[i][j] else 1 for j in range(len(prediction[i]))]

        def computeEfficiency(difference):
            eff = 0
            for predict in difference:
                eff += 1
                for feature in predict:
                    if feature != 0:
                        eff -= 1
                        break
                    
            eff /= len(difference)
            return eff*100

        self.x = self.data[:,:-4]
        self.y = self.data[:,-4:]
        pred = self.model.predict(self.x)
        process(pred)
        diff = abs(pred-self.y)
        self.eff = computeEfficiency(diff)
        print(self.name," efficiency :",self.eff,"%")
