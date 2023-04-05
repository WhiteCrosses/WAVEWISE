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
from PyQt5.QtCore import Qt
import csv

import pyqtgraph as pg

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
        self.histMin = -1
        self.histMax = 10
        
        self.colorMeshMin = self.histMin
        self.colorMeshMax = self.histMax
        
        self.bgColor = '#0b213b'
        self.setFixedWidth(1024)
        self.setFixedHeight(768)
        #self.setStyleSheet("background-color: #0b213b;")
        
        self.plutoInit()
        self.peakFinderInit()
        self.mainWidgetsInit()
        
        self.setCentralWidget(self._main)
        
        self.createWidgets()
        self.firstIteration()
        
    def firstIteration(self):
        self.x = np.arange(self.startFreq,self.bufferSize+self.startFreq)
        self.signal = self.x * 0
        self.imageArray = np.zeros((100,self.bufferSize))#change to variables later
        
        self._line = self._waterfall_ax.pcolorfast(np.reshape(self.imageArray,(100,self.bufferSize)),vmin=self.colorMeshMin,vmax=self.colorMeshMax)
        self._waterfall_ax.invert_yaxis()

    def mainWidgetsInit(self):
        self._main = QtWidgets.QWidget()
        self._leftBox = QtWidgets.QWidget()
        self._middleBox = QtWidgets.QWidget()

    def peakFinderInit(self):
        self.pSetHeight = None
        self.pSetThereshold = None
        self.pSetDistance = None
        self.pSetProminence = None
        self.pSetWidth = 1
        self.pSetWlen = 10000
        self.pSetRelHeight = 1
        self.pSetPlateauSize = None

    def plutoInit(self):
        self.constantPart = 0
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
        
    def _update_canvas(self):
        self.getData()
        self._wave_ax.clear()
        self._wave_ax.set_ylim(5,10)
        self._wave_ax.plot(self.freq,self.data)
        
        self.renderPeaks()
            
        self._line.set_array(np.reshape(self.imageArray,(100,self.bufferSize)))
        self._line.set(clim=(self.colorMeshMin,self.colorMeshMax))

        self._wave_ax.figure.canvas.draw()
        self._line.figure.canvas.draw()
    
    def renderPeaks(self):
        x,y = self.checkPeaks()
        for i in x:
            self._wave_ax.axvline(self.freq[i], color='r')
    
    def getFrequencyArray(self):
        self.freq = np.fft.fftfreq(self.signal.size,d=1/self.sampleRate)
        self.freq = self.freq + self.center_freq
        
    def getData(self):
        self.signal = self.sdr.rx()
        self.processData()
        self.getFrequencyArray()

    def dataFilter(self):
        kernel = np.ones(10)/10
        self.data = np.convolve(self.data,kernel,mode='same')
    
    def processData(self):
        
        #process pulled data from pluto
        #Fast Fourier transform
        self.data = np.fft.fft(self.signal)
        self.data = np.abs(self.data)
        
        #Converting to logarythmic scale
        self.data = np.log10(self.data*1000)
        
        #Smooth out noise
        self.dataFilter()
        
        self.data = self.data + self.constantPart
        self.addToImageArray()
        
    def addToImageArray(self):    
        self.imageArray = np.roll(self.imageArray,1,axis=0)
        self.imageArray[0, :] = np.roll(self.data,int(self.data.size/2),axis=0)
        
    def checkPeaks(self):
        x,y = signal.find_peaks(
            self.data,
            height=self.pSetHeight, 
            threshold=self.pSetThereshold, 
            distance=self.pSetDistance, 
            #prominence=self.pSetProminence, 
            width=self.pSetWidth, 
            #wlen=self.pSetWlen, 
            rel_height=self.pSetRelHeight,
            plateau_size=self.pSetPlateauSize)
        return x,y
    
    def setSampleRate(self):
        value = self.sampleRateBox.value()
        if value < 1:
            value = 1
            self.sampleRateBox.value = 1
        value = value * 1e6
        if value <= 0:
            pass
        self.sdr.sample_rate = self.sampleRate
        self.getFrequencyArray()
    
    def setCenterFreq(self):
        value = float(self.centerFreqBox.displayText())
        if value == 0:
            pass
        else:
            self.center_freq = int(value * 1e6)
            self.sdr.rx_lo = self.center_freq
        
        self.getFrequencyArray()
        
    def changeHistLevels(self):
        self.colorMeshMin, self.colorMeshMax = self.hist.getLevels()

    def pSetHeightChange(self):
        self.pSetHeight = self.pSetHeightSlider.value()
        print(self.pSetHeight)
        
    def pSetDistanceChange(self):
        self.pSetDistance = int(self.pSetDistanceSlider.value()*self.bufferSize/100)
        if self.pSetDistance < 1:
            self.pSetDistance = 1
        print(self.pSetDistance)
        
    def pSetTheresholdChange(self):
        self.pSetThereshold = self.pSetTheresholdSlider.value()/100
        print(self.pSetThereshold)
    
    def pSetProminenceChange(self):
        self.pSetProminence = self.pSetProminenceSlider.value()
        
    def pSetWidthChange(self):
        self.pSetWidth = int(self.sampleRate * self.pSetWidthSlider.value() / 10000 / 2)
        print(self.pSetWidth)
        
    def pSetWlenChange(self):
        self.pSetWlen = self.pSetWlenSlider.value()
        
    def pSetRelHeightChange(self):
        self.pSetRelHeight = self.pSetRelHeightSlider.value()
        
    def pSetPlateauSizeChange(self):
        self.pSetPlateauSize = self.pSetPlateauSizeSlider.value()

    def createWidgets(self):
        widget = QtWidgets.QWidget()
        
        button1 = QtWidgets.QPushButton(widget)
        button1.setText("Change center freq!")
        
        
        self.sampleRateBox = QtWidgets.QLineEdit()
        self.sampleRateBox.setInputMask("999.9999")
        self.sampleRateBox.setText("10.0000") 
        self.sampleRateBox.returnPressed.connect(self.setSampleRate)
        
        self.centerFreqBox = QtWidgets.QLineEdit()
        self.centerFreqBox.setInputMask("999.9999")
        self.centerFreqBox.setText("100.0000") 
        self.centerFreqBox.returnPressed.connect(self.setCenterFreq)
        
        self.hist = pg.HistogramLUTWidget(gradientPosition="left")
        self.hist.setLevels(min=self.histMin, max=self.histMax)
        self.hist.sigLevelsChanged.connect(self.changeHistLevels)
        self.hist.disableAutoHistogramRange()
        
        #Peak settings
        #height
        self.pSetHeightSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetHeightSlider.setFocusPolicy(Qt.StrongFocus)
        self.pSetHeightSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetHeightSlider.setTickInterval(1)
        self.pSetHeightSlider.setMinimum(1)
        self.pSetHeightSlider.setMaximum(10)
        self.pSetHeightSlider.setSingleStep(1)
        self.pSetHeightSlider.valueChanged.connect(self.pSetHeightChange)
        self.pSetHeightSlider.setToolTip("Required height of peaks.")
        
        #thereshold
        self.pSetTheresholdSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetTheresholdSlider.setFocusPolicy(Qt.StrongFocus)
        self.pSetTheresholdSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetTheresholdSlider.setTickInterval(1)
        self.pSetTheresholdSlider.setMinimum(1)
        self.pSetTheresholdSlider.setMaximum(10)
        self.pSetTheresholdSlider.setSingleStep(1)
        self.pSetTheresholdSlider.valueChanged.connect(self.pSetTheresholdChange)
        self.pSetTheresholdSlider.setToolTip("Required threshold of peaks, the vertical distance to its neighboring samples.")
        
        #distance
        self.pSetDistanceSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetDistanceSlider.setFocusPolicy(Qt.StrongFocus)
        self.pSetDistanceSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetDistanceSlider.setTickInterval(10)
        self.pSetDistanceSlider.setMinimum(1)
        self.pSetDistanceSlider.setMaximum(100)
        self.pSetDistanceSlider.setSingleStep(1)
        self.pSetDistanceSlider.valueChanged.connect(self.pSetDistanceChange)
        self.pSetDistanceSlider.setToolTip("Required minimal horizontal distance (>= 1) in samples between neighbouring peaks.\nSmaller peaks are removed first until the condition is fulfilled for all remaining peaks.")
        
        #prominence
        self.pSetProminenceSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetProminenceSlider.setFocusPolicy(Qt.StrongFocus)
        self.pSetProminenceSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetProminenceSlider.setTickInterval(10)
        self.pSetProminenceSlider.setSingleStep(1)
        self.pSetProminenceSlider.valueChanged.connect(self.pSetProminenceChange)
        
        #width
        self.pSetWidthSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetWidthSlider.setFocusPolicy(Qt.StrongFocus)
        self.pSetWidthSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetWidthSlider.setTickInterval(10)
        self.pSetWidthSlider.setMinimum(1)
        self.pSetWidthSlider.setMaximum(100)
        self.pSetWidthSlider.setSingleStep(1)
        self.pSetWidthSlider.valueChanged.connect(self.pSetWidthChange)
        self.pSetWidthSlider.setToolTip("Required width of peaks in samples.")
        
        #wlen
        self.pSetWlenSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetWlenSlider.setFocusPolicy(Qt.StrongFocus)
        self.pSetWlenSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetWlenSlider.setTickInterval(10)
        self.pSetWlenSlider.setSingleStep(1)
        self.pSetWlenSlider.valueChanged.connect(self.pSetWlenChange)
        
        #rel height
        self.pSetRelHeightSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetRelHeightSlider.setFocusPolicy(Qt.StrongFocus)
        self.pSetRelHeightSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetRelHeightSlider.setTickInterval(10)
        self.pSetRelHeightSlider.setSingleStep(1)
        self.pSetRelHeightSlider.valueChanged.connect(self.pSetRelHeightChange)
        self.pSetRelHeightSlider.setToolTip("Used for calculation of the peaks width, thus it is only used if width is given.")
        
        
        #plateau size
        self.pSetPlateauSizeSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetPlateauSizeSlider.setFocusPolicy(Qt.StrongFocus)
        self.pSetPlateauSizeSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetPlateauSizeSlider.setTickInterval(10)
        self.pSetPlateauSizeSlider.setSingleStep(1)
        self.pSetPlateauSizeSlider.valueChanged.connect(self.pSetPlateauSizeChange)
        self.pSetPlateauSizeSlider.setToolTip("Required size of the flat top of peaks in samples.")
        
        
        layoutpSetBox = QtWidgets.QFormLayout()
        layoutpSetBox.addRow("Height", self.pSetHeightSlider)
        layoutpSetBox.addRow("Thereshold", self.pSetTheresholdSlider)
        layoutpSetBox.addRow("Distance", self.pSetDistanceSlider)
        #layoutpSetBox.addRow("Prominence", self.pSetProminenceSlider)
        layoutpSetBox.addRow("Width", self.pSetWidthSlider)
        #layoutpSetBox.addRow("Wlen", self.pSetWlenSlider)
        layoutpSetBox.addRow("Relative Height", self.pSetRelHeightSlider)
        layoutpSetBox.addRow("Plateau Size", self.pSetPlateauSizeSlider)


        pSetBoxWidget = QtWidgets.QGroupBox("Peak settings")
        pSetBoxWidget.setLayout(layoutpSetBox)
        
        layoutMainBox = QtWidgets.QHBoxLayout(self._main)
        layoutMiddleBox = QtWidgets.QVBoxLayout()
        layoutLeftBox = QtWidgets.QVBoxLayout()
        
        
        wave_canvas = FigureCanvas(Figure(figsize=(10, 6)))
        waterfall_canvas = FigureCanvas(Figure(figsize=(10, 6)))
        
        
        layoutMiddleBox.addWidget(wave_canvas)
        layoutMiddleBox.addWidget(waterfall_canvas)
        
        
        layoutSamplingSetBox = QtWidgets.QFormLayout()
        layoutSamplingSetBox.addRow("Sample rate", self.sampleRateBox)
        layoutSamplingSetBox.addRow("Center Frequency", self.centerFreqBox)
        samplingSetBoxWidget = QtWidgets.QGroupBox("Peak settings")
        samplingSetBoxWidget.setLayout(layoutSamplingSetBox)
        
        
        layoutLeftBox.addWidget(samplingSetBoxWidget)
        layoutLeftBox.addWidget(pSetBoxWidget)
        
        
        self._leftBox.setLayout(layoutLeftBox)
        self._middleBox.setLayout(layoutMiddleBox)
        
        
        layoutMainBox.addWidget(self._leftBox)
        layoutMainBox.addWidget(self._middleBox)
        layoutMainBox.addWidget(self.hist)
        
        #button1.clicked.connect(self.moveCenter)
        
        wave_canvas.figure.set_facecolor(self.bgColor)
        waterfall_canvas.figure.set_facecolor(self.bgColor)
        
        with plt.rc_context({'axes.edgecolor':'white', 'xtick.color':'white', 'ytick.color':'white', 'figure.facecolor':self.bgColor}):
            self._wave_ax = wave_canvas.figure.subplots()
            self._waterfall_ax = waterfall_canvas.figure.subplots()
            self._wave_ax.tick_params(colors='white',which='both')
            
        self._timer = waterfall_canvas.new_timer(100)
        self._timer.add_callback(self._update_canvas)
        self._timer.start()
    
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