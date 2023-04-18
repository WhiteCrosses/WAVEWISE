import sys
import time
import adi

import numpy as np



from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import pandas as pd
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QMainWindow, QAction, qApp, QApplication

from PyQt5.QtGui import QIcon
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
        self.histMin = -100
        self.histMax = 0
        
        self.constantPart = 0
        self.sampleRate = int(10e6)
        self.center_freq = int(100e6)
        self.bufferSize = 1024
        self.startFreq = self.center_freq-(self.sampleRate/2)
        
        self.frameTime = 100
        
        self.isRecording = False
        self.isFirstIteration = False
        self.isPlutoRunning = False
        self.isFile = False
        
        self.filterFrame = 5
        self.filterFrameLength = 20
        
        self.colorMeshMin = self.histMin
        self.colorMeshMax = self.histMax
        
        self.bgColor = '#0b213b'
        self.setFixedWidth(1600)
        self.setFixedHeight(900)
        
        self.peaks = []
        
        #self.setStyleSheet("background-color: #0b213b;")
        
       
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

    def selectFle(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', 
         'c:\\',"Text files (*.csv *.txt)")
        
        print(f"opened: {fname}")


    def peakFinderInit(self):
        self.pSetHeight = -50
        self.pSetThereshold = None
        self.pSetDistance = 10
        self.pSetProminence = None
        self.pSetWidth = 1
        self.pSetWlen = 10000
        self.pSetRelHeight = 1
        self.pSetPlateauSize = None

    def plutoInit(self):
        self.sdr = adi.Pluto("ip:192.168.2.1")
        self.sdr.sample_rate = self.sampleRate
        self.sdr.rx_rf_bandwidth = self.sampleRate # filter cutoff, just set it to the same as sample rate
        self.sdr.rx_lo = self.center_freq
        self.sdr.rx_buffer_size = self.bufferSize # this is the buffer the Pluto uses to buffer samples
        self.sdr.rx_rf_bandwidth=int(self.sampleRate)
        
    def _update_canvas(self):
        
        if self.isFirstIteration:
             self.plutoInit()
             self.isFirstIteration = False
        
        if self.isPlutoRunning:
            self.getData()
            self._wave_ax.clear()
            self._wave_ax.set_ylim(-100,0)
        
            self._wave_ax.set_ylabel("dBm", color='white')
        
            def millions(x, pos):
                if x<1e6:
                    return '%1.1fkHz' % (x*1e-3)
                elif x>=1e6 and x<1e9:
                    return '%1.1fMHz' % (x*1e-6)
                elif x>=1e9:
                    return '%1.1fGHz' % (x*1e-9)
                

            formatter = FuncFormatter(millions)
            self._wave_ax.xaxis.set_major_formatter(formatter)
        
            self._wave_ax.plot(self.freq,self.data)
            
            self.renderPeaks()
                
            self._line.set_array(np.reshape(self.imageArray,(100,self.bufferSize)))
            self._line.set(clim=(self.colorMeshMin,self.colorMeshMax))

            self._wave_ax.figure.canvas.draw()
            self._line.figure.canvas.draw()
        
            self.peaks.sort(key=lambda x:x[0])
            for row in range(len(self.peaks)):
                self.pRapTable.setItem(row,0,QtWidgets.QTableWidgetItem(str(self.peaks[row][0])))
                print(f"Adding {self.peaks[row][0]} at: {row}")
        
    def renderPeaks(self):
        x,y = self.checkPeaks()
        buffer_size = self.pSetDistance / 2
        print(buffer_size)
        peak_array = x
        
        
        
        for n in peak_array:
            peakExists = False
            for m in self.peaks:
                #print(f"n = {n}\nm = {m[0]}")
                m[1]=False
                if n+buffer_size > m[0] > n-buffer_size:
                    peakExists = True
                    
                    print(f"existing peak found at: {n}, delta: {np.abs(m[0]-n)}")
                    m = [n,True]
                    break
                #print("\n")
            
            if not peakExists:
                self.peaks.append([n,True,3])

        for m in self.peaks:
            if m[1]==False and m[2]<=0:
                print(f"removed:{m}")
                self.peaks.remove(m)
            elif m[1]==False and m[2]>0:
                m[2]-=1
            elif m[1]==True:
                m[2]=5
        
        print("Peaks:")
        print(self.peaks)
                     
        for i in x:
            self._wave_ax.axvline(self.freq[i], color='r')
    
    def getFrequencyArray(self):
        self.freq = np.fft.fftfreq(self.signal.size,d=1/self.sampleRate)
        self.freq = self.freq + self.center_freq
        
    def getData(self):
        if self.isPlutoRunning:
            self.signal = self.sdr.rx()
            self.processData()
            self.getFrequencyArray()

    def dataFilter(self):
        kernel = np.ones(10)/10
        self.data = np.convolve(self.data,kernel,mode='same')
    
    def processData(self):
        
        #process pulled data from pluto
        #Fast Fourier transform
        self.data = np.abs(np.fft.fft(self.signal))**2 / (self.bufferSize*self.sampleRate)
        
        #Converting to logarythmic scale
        self.data = (20 * np.log10(np.abs(self.data)))
        
        thereshold = np.average(self.data)*0.9
        for i in range(len(self.data)):
            if self.data[i] < thereshold:
                self.data[i] = np.average(self.data)
        #Smooth out noise
        self.dataFilter()
        
        #self.data = self.data + self.constantPart
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
        
    def plutoStartChangeState(self):
        self.isPlutoRunning = not self.isPlutoRunning
        self.isFirstIteration = not self.isFirstIteration
        
    def changeHistLevels(self):
        self.colorMeshMin, self.colorMeshMax = self.hist.getLevels()

    """
    Methods reacting to changes in peak finder settings.
    """
    def pSetHeightChange(self):
        self.pSetHeight = self.pSetHeightSlider.value()
        self.pSetHeightSpinBox.setValue(self.pSetHeight)
        print(self.pSetHeight)
  
    def pSetHeightChangeSpinBox(self):
        self.pSetHeight = self.pSetHeightSpinBox.value()
        self.pSetHeightSlider.setValue(self.pSetHeight)
        print(self.pSetHeight)
              
    def pSetDistanceChange(self):
        self.pSetDistance = int(self.pSetDistanceSlider.value()*self.bufferSize/100)
        if self.pSetDistance < 1:
            self.pSetDistance = 1
        print(self.pSetDistance)
        
    def pSetTheresholdChange(self):
        self.pSetThereshold = self.pSetTheresholdSlider.value() / 100
        self.pSetTheresholdSpinBox.setValue(int(self.pSetThereshold * 10))
        print(self.pSetThereshold)
    
    def pSetTheresholdChangeSpinBox(self):
        self.pSetThereshold = self.pSetTheresholdSpinBox.value() / 100
        self.pSetTheresholdSlider.setValue(int(self.pSetThereshold * 10))
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

    def recChangeRecordingState(self):
        
        if self.isFile:
            self.isRecording = not self.isRecording
            self.file
            if self.isRecording:
                self.recButton.setText("Stop")
            elif not self.isRecording:
                self.recButton.setText("Start")
        else:
            self.selectFle()
        
    def createWidgets(self):
        widget = QtWidgets.QWidget()
        self.menuBar = QtWidgets.QMenuBar(self)
        self.setMenuBar(self.menuBar)
        
        
        
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
        self.pSetHeightSlider.setTickInterval(5)
        self.pSetHeightSlider.setMinimum(-100)
        self.pSetHeightSlider.setMaximum(0)
        self.pSetHeightSlider.setSingleStep(5)
        self.pSetHeightSlider.valueChanged.connect(self.pSetHeightChange)
        self.pSetHeightSlider.setToolTip("Required height of peaks.")
        
        
        self.pSetHeightSpinBox = QtWidgets.QSpinBox()
        self.pSetHeightSpinBox.valueChanged.connect(self.pSetHeightChangeSpinBox)
        self.pSetHeightSpinBox.setRange(-100,0)
        self.pSetHeightSpinBox.setSingleStep(5)
        
        
        self.pSetHeightSlider.setValue(-50)
        self.pSetHeightSpinBox.setValue(-50)
        
        
        self.pSetHeightBoxLayout = QtWidgets.QHBoxLayout()
        self.pSetHeightBoxLayout.addWidget(self.pSetHeightSlider)
        self.pSetHeightBoxLayout.addWidget(self.pSetHeightSpinBox)
        self.pSetHeightBox = QtWidgets.QHBoxLayout()
        self.pSetHeightBox.addLayout(self.pSetHeightBoxLayout)
        
        
        #thereshold
        self.pSetTheresholdSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetTheresholdSlider.setFocusPolicy(Qt.StrongFocus)
        self.pSetTheresholdSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetTheresholdSlider.setTickInterval(5)
        self.pSetTheresholdSlider.setMinimum(1)
        self.pSetTheresholdSlider.setMaximum(15)
        self.pSetTheresholdSlider.setSingleStep(1)
        self.pSetTheresholdSlider.valueChanged.connect(self.pSetTheresholdChange)
        self.pSetTheresholdSlider.setToolTip("Required threshold of peaks, the vertical distance to its neighboring samples.")
        
        
        self.pSetTheresholdSpinBox = QtWidgets.QSpinBox()
        self.pSetTheresholdSpinBox.valueChanged.connect(self.pSetTheresholdChangeSpinBox)
        self.pSetTheresholdSpinBox.setRange(1,15)
        self.pSetTheresholdSpinBox.setSingleStep(1)
        
        
        self.pSetTheresholdSlider.setValue(2)
        self.pSetTheresholdSpinBox.setValue(2)
        
        
        self.pSetTheresholdBoxLayout = QtWidgets.QHBoxLayout()
        self.pSetTheresholdBoxLayout.addWidget(self.pSetTheresholdSlider)
        self.pSetTheresholdBoxLayout.addWidget(self.pSetTheresholdSpinBox)
        self.pSetTheresholdBox = QtWidgets.QHBoxLayout()
        self.pSetTheresholdBox.addLayout(self.pSetTheresholdBoxLayout)
        
        
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
        
        
        #Peak tracking table
        self.pRapTable = QtWidgets.QTableWidget(10,3)
        
        
        #Start/Stop button for recording to csv
        self.recButton = QtWidgets.QPushButton("Start")
        self.recButton.clicked.connect(self.recChangeRecordingState)
        recButtonLayout = QtWidgets.QHBoxLayout()
        recButtonLayout.addWidget(self.recButton)
        recordingBoxWidget = QtWidgets.QGroupBox("Recording")
        recordingBoxWidget.setLayout(recButtonLayout)
        
        
        self.plutoStartButton = QtWidgets.QPushButton("Start Pluto!")
        self.plutoStartButton.clicked.connect(self.plutoStartChangeState)
        plutoStartButtonLayout = QtWidgets.QHBoxLayout()
        plutoStartButtonLayout.addWidget(self.plutoStartButton)
        plutoStartBoxWidget = QtWidgets.QGroupBox("I/O")
        plutoStartBoxWidget.setLayout(plutoStartButtonLayout)
        
        
        layoutpSetBox = QtWidgets.QFormLayout()
        layoutpSetBox.addRow("Height", self.pSetHeightBox)
        layoutpSetBox.addRow("Thereshold", self.pSetTheresholdBox)
        layoutpSetBox.addRow("Distance", self.pSetDistanceSlider)
        #layoutpSetBox.addRow("Prominence", self.pSetProminenceSlider)
        layoutpSetBox.addRow("Width", self.pSetWidthSlider)
        #layoutpSetBox.addRow("Wlen", self.pSetWlenSlider)
        layoutpSetBox.addRow("Relative Height", self.pSetRelHeightSlider)
        layoutpSetBox.addRow("Plateau Size", self.pSetPlateauSizeSlider)

        layoutpRapBox = QtWidgets.QTableWidget

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
        
        
        
        
        
        layoutLeftBox.addWidget(plutoStartBoxWidget)
        layoutLeftBox.addWidget(recordingBoxWidget)
        layoutLeftBox.addWidget(samplingSetBoxWidget)
        layoutLeftBox.addWidget(pSetBoxWidget)
        layoutLeftBox.addWidget(self.pRapTable)
        
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