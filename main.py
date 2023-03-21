import sys
import time
import adi

import numpy as np

from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import pandas as pd
import csv

from scipy import signal
from scipy.fft import fftshift


class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        
        """
        
            Signal processing setting variables:
                freqRange - int, Hz
                startFreq - int, Hz
                sampleRate - int, no. samples per second
                selfDuration - int, seconds, time of collecting samples
        
        """

        self.selectedMode = 1
        
        
        self.sineArray = np.array([[1,828]])#[amplitude,frequency]
        widget = QtWidgets.QWidget()
        
        self.deviceSelectionBox = QtWidgets.QComboBox(self)
        self.deviceSelectionBox.addItem("Adalm-Pluto")
        self.deviceSelectionBox.addItem("Emulator")
        
        
        button1 = QtWidgets.QPushButton(widget)
        button1.setText("Change center freq!")
        
        self.sineAmpBox = QtWidgets.QLineEdit(self)
        self.sineAmpBox.setText("1351")
        
        self.x = np.arange(self.startFreq,self.bufferSize+self.startFreq)
        print(self.x)
        self.signal = self.x * 0
        self.imageArray = np.zeros((100,self.bufferSize))#change to variables later
        
        self._main = QtWidgets.QWidget()
        self._leftBox = QtWidgets.QWidget()
        self._middleBox = QtWidgets.QWidget()
        
        self.setCentralWidget(self._main)
        
        layoutMainBox = QtWidgets.QHBoxLayout(self._main)
        layoutMiddleBox = QtWidgets.QVBoxLayout()
        layoutLeftBox = QtWidgets.QVBoxLayout()
        
        
        wave_canvas = FigureCanvas(Figure(figsize=(5, 3)))
        waterfall_canvas = FigureCanvas(Figure(figsize=(5, 3)))
        
        layoutMiddleBox.addWidget(wave_canvas)
        layoutMiddleBox.addWidget(waterfall_canvas)
        
        layoutLeftBox.addWidget(self.deviceSelectionBox)
        layoutLeftBox.addWidget(button1)
        layoutLeftBox.addWidget(self.sineAmpBox)
        
        
        self._leftBox.setLayout(layoutLeftBox)
        self._middleBox.setLayout(layoutMiddleBox)
        
        layoutMainBox.addWidget(self._leftBox)
        layoutMainBox.addWidget(self._middleBox)
        
        button1.clicked.connect(self.moveCenter)
        self._wave_ax = wave_canvas.figure.subplots()
        self._waterfall_ax = waterfall_canvas.figure.subplots()

        self.gnerateData()
        self._line = self._waterfall_ax.pcolorfast(np.reshape(self.imageArray,(100,self.bufferSize)),vmin=np.min(self.imageArray),vmax=np.max(self.imageArray))
        self._waterfall_ax.invert_yaxis()
        #self._wave_ax.plot(0,0)
        self._wave_ax.autoscale(False)
        
        self._timer = waterfall_canvas.new_timer(100)
        self._timer.add_callback(self._update_canvas)
        self._timer.start()

    def plutoInit(self):
        self.sampleRate = int(10e6)
        self.center_freq = int(100e6)
        self.bufferSize = 1024
        self.startFreq = self.center_freq-(self.sampleRate/2)
        
        self.sdr = adi.Pluto("ip:192.168.2.1")
        self.sdr.sample_rate = self.sampleRate
        self.sdr.rx_rf_bandwidth = self.sampleRate # filter cutoff, just set it to the same as sample rate
        self.sdr.rx_lo = self.center_freq
        self.sdr.rx_buffer_size = self.bufferSize # this is the buffer the Pluto uses to buffer samples
        self.sdr.rx_rf_bandwidth=int(self.sampleRate)

    def moveCenter(self):
        self.center_freq = self.center_freq + int(5e6)
        self.sdr.rx_lo = self.center_freq

    def _update_canvas(self):
        
        ampVal = self.sineAmpBox.text()
        self.selectedMode = self.deviceSelectionBox.currentIndex()
        try:
            ampVal = int(ampVal)
            self.sineArray[0][1] = ampVal
        except ValueError:
            pass
        self.getData()
        self._wave_ax.clear()
        self._wave_ax.plot(self.freq,self.data)
        
        x,y = self.checkPeaks()
        #for i in x:
            #self._wave_ax.axvline(x=i, color='r')
            
        self._line.set_array(np.reshape(self.imageArray,(100,self.bufferSize)))
        self._line.set(clim=(np.min(self.imageArray),np.max(self.imageArray)))

        self._wave_ax.figure.canvas.draw()
        self._line.figure.canvas.draw()
        print(np.min(self.freq))
        print(np.max(self.freq))
    
        
    def getData(self):
        if self.selectedMode == 0:
            self.signal = self.sdr.rx()
        else:
            self.gnerateData()
        self.processData()
        
        self.freq = np.fft.fftfreq(self.signal.size,d=1/self.sampleRate)
        self.freq = self.freq + self.center_freq
        print(self.freq)

    def shortenSignal(self):
        averages = [sum(self.data[i:i+5])/5 for i in range (0, len(self.data),5)]
        self.data = averages
        
    def gnerateData(self):
        self.freqArray = np.fft.rfftfreq(self.sampleRate, 1/self.sampleRate)
        self.signal = np.sum(np.apply_along_axis(self.createSine, axis=1,arr=self.sineArray),axis=0)
        
        
    def processData(self):
    
        self.data = np.fft.fft(self.signal)
        self.data = np.abs(self.data)
        #self.freqArray = self.freqArray[:-1]
        self.addToImageArray()
        
        
    def addToImageArray(self):    
        self.imageArray = np.roll(self.imageArray,1,axis=0)
        self.imageArray[0, :] = self.data

        
    def createSine(self,array):
        x = array[0] * np.sin(2*np.pi*self.x*array[1])
        return x 
        
        
    def checkPeaks(self):
        x,y = signal.find_peaks(self.data,height=(np.max(self.data)/3))
        return x,y
    
    
    def addSine(self):
        newSine = self.sineArray[0:]
        delta = [[1,5]]
        newSine = np.add(newSine,delta)
        
        self.sineArray = np.concatenate((self.sineArray,newSine))
        
        

if __name__ == "__main__":
    # Check whether there is already a running QApplication (e.g., if running
    # from an IDE).
    qapp = QtWidgets.QApplication.instance()
    if not qapp:
        qapp = QtWidgets.QApplication(sys.argv)

    app = ApplicationWindow()
    app.show()
    app.activateWindow()
    app.raise_()
    qapp.exec()