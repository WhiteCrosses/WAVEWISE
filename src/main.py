import pyqtgraph as pg

import sys
import adi

import numpy as np

#Krzywe transmisyjne

from superqt import QLabeledRangeSlider, QLabeledSlider

import matplotlib.transforms as transforms
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt, QEvent

import csv



from scipy import signal

class Peak:
    def __init__(self, frequency, power, distance):
        self.frequency = frequency
        self.distance = distance
        self.min = self.frequency - self.distance
        self.max = self.frequency + self.distance
        self.tickCounter = 0
        self.isChecked = False
        self.power = power

    def found(self, x, power):
        self.frequency = x
        self.power = power
        self.min = self.frequency - self.distance
        self.max = self.frequency + self.distance
        self.resetTickCounter()
        
    def resetTickCounter(self):
        self.isChecked = True
        self.tickCounter = 0

#foo bar
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
        
        
        self.recType = 0
        self.peakArray = []
        self.delay = 100
        
        self.constantPart = 0
        self.sampleRate = int(10e6)
        self.center_freq = int(100e6)
        self.bufferSize = 1024
        self.startFreq = self.center_freq-(self.sampleRate/2)
        
        self.recMarkerValues = (0,0)
        
        self.isPlutoRunning = False
        self.frameTime = 100
        
        self.isRecording = False
        self.isFirstIteration = False
        
        self.isFile = False
        
        self.filterFrame = 5
        self.filterFrameLength = 20
        
        self.colorMeshMin = self.histMin
        self.colorMeshMax = self.histMax
        
        self.bgColor = '#0b213b'
        self.setGeometry(0,0,1600,900)
        
        #self.setFixedWidth(1600)
        #self.setFixedHeight(900)
        
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
        ticks = np.array([0,1,2,3,4])
        _ticks = ticks/4*self.sampleRate+self.startFreq
        
        
        self.calculateWaterfallNodes()
        #self._waterfall_ax.xaxis.set_major_locator(ticker.FixedLocator(_ticks))
        
        self._wave_ax.set_ylim(-100,0)
        self._wave_ax.set_ylim(0,100)
    
    def mainWidgetsInit(self):
        self._main = QtWidgets.QWidget()
        self._leftBox = QtWidgets.QWidget()
        self._middleBox = QtWidgets.QWidget()

    def selectFle(self):
        
        print("opened:")

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
                    return '%1.4fGHz' % (x*1e-9)
                

            formatter = FuncFormatter(millions)
            self._wave_ax.xaxis.set_major_formatter(formatter)
            #self.freq = np.sort(self.freq)
            #print(self.freq)
            
            
            dataarray = list(zip(self.freq, self.data))
            dataarray = sorted(dataarray, key=lambda x: x[0])
            
            self.freq, self.data = zip(*dataarray)
            
            self._wave_ax.plot(self.freq, self.data, color='blue')
            
            self.renderPeaks()
            
            self._line.set_array(np.reshape(self.imageArray,(100,self.bufferSize)))
            self._line.set(clim=(self.colorMeshMin,self.colorMeshMax))

            self._wave_ax.figure.canvas.draw()
            self._line.figure.canvas.draw()
        
            for row in range(10):
                self.pRapTable.setItem(row,0,QtWidgets.QTableWidgetItem(""))
            
            for row in range(len(self.peakArray)):
                self.pRapTable.setItem(row,0,QtWidgets.QTableWidgetItem(str(self.freq[self.peakArray[row].frequency])))
                 
    def renderPeaks(self):
        
        x,y = self.checkPeaks()  
        
        #Iterate through peaks to assign them new positions
        m=0
        for i in x:
            found = False
            for j in self.peakArray:
                if j.min<=i<=j.max:
                    j.found(i,y['peak_heights'][m])
                    found = True
                    break
                
            if not found:
                self.peakArray.append(Peak(i,y['peak_heights'][m],self.pSetDistance))
            m+=1
        
        for j in self.peakArray:
            if not j.isChecked:
                j.tickCounter += 1
                if j.tickCounter > 10:
                    self.peakArray.remove(j)
        
        for j in self.peakArray:
            j.isChecked = False
                     
        for i in x:
            self._wave_ax.axvline(self.freq[i], color='r')
              
    def getFrequencyArray(self):
        self.freq = np.fft.fftfreq(1024,d=1/self.sampleRate)
        self.freq = self.freq + self.center_freq
        
    def getData(self):
        if self.isPlutoRunning:
            try:
                self.signal = self.sdr.rx()
                self.processData()
                self.getFrequencyArray()
            except:
                pass

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
    
    def calculateWaterfallNodes(self):
        
        def millions(x):
                if x<1e6:
                    return '%1.1fkHz' % (x*1e-3)
                elif x>=1e6 and x<1e9:
                    return '%1.1fMHz' % (x*1e-6)
                elif x>=1e9:
                    return '%1.4fGHz' % (x*1e-9)
        
        
        
        ticks = np.array([0,1,2,3,4])
        _ticks = ticks/4*self.sampleRate+self.startFreq
        
        formattedTicks = ['','','','','']
        for i in range(len(_ticks)):
            formattedTicks[i] = millions(_ticks[i])
        
        tickLoc = np.arange(0,1024,1024/4)
        tickLoc = np.append(tickLoc,1024)
        print(tickLoc)
        self._waterfall_ax.set_xticks(tickLoc)
        
        self._waterfall_ax.locator_params(axis='both', nbins=5)
        self._waterfall_ax.set_xticklabels(formattedTicks)
    
    def setSampleRate(self):
        value = self.sampleRateBox.text()
        value = float(value)
        if value < 1:
            value = 1
            self.sampleRateBox.value = 1
        value = value * 1e6
        if value <= 0:
            pass
        self.sampleRate = value
        self.sdr.sample_rate = self.sampleRate
        self.plutoInit()
        self.getFrequencyArray()

        print("sample rate has been changed")
        
    def setCenterFreq(self):
        value = float(self.centerFreqBox.displayText())
        print(value)
        if value <= 0:
            pass
        else:
            self.center_freq = int(value * 1e6)
            self.startFreq = self.center_freq-(self.sampleRate/2)
            self.sdr.rx_lo = self.center_freq
        
        self.calculateWaterfallNodes()
        
        
        self.getFrequencyArray()
        
    def plutoStartChangeState(self):
        self.isPlutoRunning = not self.isPlutoRunning
        if self.isPlutoRunning:
            self.plutoStartButton.setText("Stop Pluto!")
        else:
            self.plutoStartButton.setText("Start Pluto!")
        self.isFirstIteration = not self.isFirstIteration
        
    def changeHistLevels(self):
        self.colorMeshMin, self.colorMeshMax = self.hist.getLevels()

    """
    Methods reacting to changes in peak finder settings.
    """
    def pSetHeightChange(self):
        self.pSetHeight = self.pSetHeightSlider.value()
        self.pSetHeightSpinBox.setValue(self.pSetHeight)
  
    def pSetHeightChangeSpinBox(self):
        self.pSetHeight = self.pSetHeightSpinBox.value()
        self.pSetHeightSlider.setValue(self.pSetHeight)
              
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
        self.recResolve()
      
    def eventFilter(self, object, event):
        if event.type() == QEvent.Enter:
            self.statusBar.showMessage(object.description)
            return True
        elif event.type() == QEvent.Leave:
            print("Mouse is not over the label")
                    
        return False
    
    def recTypeChange(self, index):
        self.recType = index
                
    def recRangeChange(self):
        self.recMarkerValues = self.recRangeSelector.value()
        
    def savePeaks(self):
        with open('profiles2.csv', 'w', newline='') as file:
            writer = csv.writer(file)            
            writer.writerow(["Frequency [Hz]", "Power [dBm]"])
            m=0
            for i in self.peakArray:
                writer.writerow([self.freq[i.frequency],i.power])
            file.close()
    
    def saveWaterfall(self):
        self._waterfall_ax.figure.savefig('waterfall.png')
    
    def recResolve(self):
        
        match self.recType:
            case 0:
                self.savePeaks()
            case 1:
                self.saveWaterfall()
            case 2:
                self.savePeaks()
                self.saveWaterfall()
    def createWidgets(self):
        widget = QtWidgets.QWidget()
        self.menuBar = QtWidgets.QMenuBar(self)
        self.setMenuBar(self.menuBar)
        
        
        
        button1 = QtWidgets.QPushButton(widget)
        button1.setText("Change center freq!")
        
        
        self.sampleRateBox = QtWidgets.QLineEdit()
        self.sampleRateBox.setInputMask("99.9999")
        self.sampleRateBox.setText("10.0000") 
        self.sampleRateBox.returnPressed.connect(self.setSampleRate)
        self.sampleRateBox.description = "Sample Rate in MHz."
        self.sampleRateBox.installEventFilter(self)
        
        self.centerFreqBox = QtWidgets.QLineEdit()
        self.centerFreqBox.setInputMask("9999.9999")
        self.centerFreqBox.setText("0100.0000") 
        self.centerFreqBox.returnPressed.connect(self.setCenterFreq)
        self.centerFreqBox.description = "Center frequency in MHz."
        self.centerFreqBox.installEventFilter(self)
        
        self.hist = pg.HistogramLUTWidget(gradientPosition="left")
        self.hist.setLevels(min=self.histMin, max=self.histMax)
        self.hist.sigLevelsChanged.connect(self.changeHistLevels)
        self.hist.disableAutoHistogramRange()
        self.hist.description = "Select range for color saturation on histogram."
        self.hist.installEventFilter(self)
        
        #Peak settings
        #height
        self.pSetHeightSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetHeightSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetHeightSlider.setTickInterval(5)
        self.pSetHeightSlider.setMinimum(-100)
        self.pSetHeightSlider.setMaximum(0)
        self.pSetHeightSlider.setSingleStep(5)
        self.pSetHeightSlider.valueChanged.connect(self.pSetHeightChange)
        #self.pSetHeightSlider.setToolTip("Required height of peaks.")
        self.pSetHeightSlider.description = "Declare minimum height to peak be detected. Values ranging from -100 to 0."
        self.pSetHeightSlider.installEventFilter(self)
        
        
        self.pSetHeightSpinBox = QtWidgets.QSpinBox()
        self.pSetHeightSpinBox.valueChanged.connect(self.pSetHeightChangeSpinBox)
        self.pSetHeightSpinBox.setRange(-100,0)
        self.pSetHeightSpinBox.setSingleStep(5)
        self.pSetHeightSpinBox.description = "Declare minimum height to peak be detected. Values ranging from -100 to 0."
        self.pSetHeightSpinBox.installEventFilter(self)
        
        
        self.pSetHeightSlider.setValue(-50)
        self.pSetHeightSpinBox.setValue(-50)
        
        
        self.pSetHeightBoxLayout = QtWidgets.QHBoxLayout()
        self.pSetHeightBoxLayout.addWidget(self.pSetHeightSlider)
        self.pSetHeightBoxLayout.addWidget(self.pSetHeightSpinBox)
        self.pSetHeightBox = QtWidgets.QHBoxLayout()
        self.pSetHeightBox.addLayout(self.pSetHeightBoxLayout)
        
        #thereshold
        self.pSetTheresholdSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetTheresholdSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetTheresholdSlider.setTickInterval(5)
        self.pSetTheresholdSlider.setMinimum(1)
        self.pSetTheresholdSlider.setMaximum(15)
        self.pSetTheresholdSlider.setSingleStep(1)
        self.pSetTheresholdSlider.valueChanged.connect(self.pSetTheresholdChange)
        #self.pSetTheresholdSlider.setToolTip("Required threshold of peaks, the vertical distance to its neighboring samples.")
        self.pSetTheresholdSlider.description = "Declare minimum diference in the values between samples"
        self.pSetTheresholdSlider.installEventFilter(self)
        
        
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
        self.pSetDistanceSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetDistanceSlider.setTickInterval(10)
        self.pSetDistanceSlider.setMinimum(1)
        self.pSetDistanceSlider.setMaximum(100)
        self.pSetDistanceSlider.setSingleStep(1)
        self.pSetDistanceSlider.valueChanged.connect(self.pSetDistanceChange)
        #self.pSetDistanceSlider.setToolTip("Required minimal horizontal distance (>= 1) in samples between neighbouring peaks.\nSmaller peaks are removed first until the condition is fulfilled for all remaining peaks.")
        self.pSetDistanceSlider.description = "Required minimal horizontal distance (>= 1) in samples between neighbouring peaks.\nSmaller peaks are removed first until the condition is fulfilled for all remaining peaks."
        self.pSetDistanceSlider.installEventFilter(self)
        
        
        #prominence
        self.pSetProminenceSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetProminenceSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetProminenceSlider.setTickInterval(10)
        self.pSetProminenceSlider.setSingleStep(1)
        self.pSetProminenceSlider.valueChanged.connect(self.pSetProminenceChange)
        
        #width
        self.pSetWidthSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetWidthSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetWidthSlider.setTickInterval(10)
        self.pSetWidthSlider.setMinimum(1)
        self.pSetWidthSlider.setMaximum(100)
        self.pSetWidthSlider.setSingleStep(1)
        self.pSetWidthSlider.valueChanged.connect(self.pSetWidthChange)
        #self.pSetWidthSlider.setToolTip(f"Required width of peaks in percentage of all samples ({self.bufferSize}).")
        self.pSetWidthSlider.description = f"Required width of peaks in percentage of all samples ({self.bufferSize})."
        self.pSetWidthSlider.installEventFilter(self)
        
        #wlen
        self.pSetWlenSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetWlenSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetWlenSlider.setTickInterval(10)
        self.pSetWlenSlider.setSingleStep(1)
        self.pSetWlenSlider.valueChanged.connect(self.pSetWlenChange)
        
        #rel height
        self.pSetRelHeightSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetRelHeightSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetRelHeightSlider.setTickInterval(10)
        self.pSetRelHeightSlider.setSingleStep(1)
        self.pSetRelHeightSlider.valueChanged.connect(self.pSetRelHeightChange)
        #self.pSetRelHeightSlider.setToolTip("Used for calculation of the peaks width, thus it is only used if width is given.")
        
        
        #plateau size
        self.pSetPlateauSizeSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetPlateauSizeSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetPlateauSizeSlider.setTickInterval(10)
        self.pSetPlateauSizeSlider.setSingleStep(1)
        self.pSetPlateauSizeSlider.valueChanged.connect(self.pSetPlateauSizeChange)
        #self.pSetPlateauSizeSlider.setToolTip("Required size of the flat top of peaks in samples.")
        
        self.statusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.statusBar)
        self.statusBar.showMessage(" ")
        
        self.recTypeChoice = QtWidgets.QComboBox()
        self.recTypeChoice.addItems(("Save peaks",
                                    "Save waterfall",
                                    "Save both"))
        self.recTypeChoice.currentIndexChanged.connect(self.recTypeChange)
        self.recTypeChoice.description = "Declare requested data to be saved"
        self.recTypeChoice.installEventFilter(self)
        
        self.recRangeSelector = QLabeledRangeSlider(Qt.Horizontal)
        self.recRangeSelector.setValue((0,121))
        self.recRangeSelector.valueChanged.connect(self.recRangeChange)
        
        #Peak tracking table
        self.pRapTable = QtWidgets.QTableWidget(10,3)
        
        
        #Start/Stop button for recording to csv
        self.recButton = QtWidgets.QPushButton("Save!")
        self.recButton.clicked.connect(self.recChangeRecordingState)
        recButtonLayout = QtWidgets.QGridLayout()
        recButtonLayout.addWidget(self.recButton,0,0)
        recButtonLayout.addWidget(self.recTypeChoice,0,1)
        self.recButton.description = "Save data to the file."
        self.recButton.installEventFilter(self)
        #   recButtonLayout.addWidget(self.recRangeSelector,1,0,1,2)
        
        recordingBoxWidget = QtWidgets.QGroupBox("Recording")
        recordingBoxWidget.setLayout(recButtonLayout)
        
        
        self.plutoStartButton = QtWidgets.QPushButton("Start Pluto!")
        self.plutoStartButton.installEventFilter(self)
        self.plutoStartButton.description = "Start the device"
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
        #layoutpSetBox.addRow("Relative Height", self.pSetRelHeightSlider)
        #layoutpSetBox.addRow("Plateau Size", self.pSetPlateauSizeSlider)

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
            
        self._timer = waterfall_canvas.new_timer(self.delay)
        self._timer.add_callback(self._update_canvas)
        self._timer.start()

    def closeEvent(self, event):
        plot.close()
        event.accept()

class TransmitWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.frequencyTable = np.array(np.zeros((3,2)))
        self.buffer = 1024
        self.gain = -50
        self.bufferAranged = np.array(np.arange(0,1024,1))
        self.mainLayout = QtWidgets.QVBoxLayout()
        self.done = True
        self._initPlotWidget()
        self._initFrequencyAddingWidget()
        
        self.setLayout(self.mainLayout)
        
    def _initPlotWidget(self):
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.mainLayout.addWidget(self.canvas)
        self.figure.clear()
        
        self.ax = self.figure.add_subplot(111)
        data = np.zeros((1024))
        
        self.ax.plot(self.bufferAranged, data, '*-')
        self.canvas.draw()

    def gainChange(self):
        self.gain = self.gainWidget.value()

    def normalize(self,arr):
        max_magnitude = np.max(np.abs(arr))  # Find the maximum absolute magnitude in the array
        normalized_arr = arr / max_magnitude  # Divide the array by the maximum magnitude
        return normalized_arr

    def runTransmit(self):
        self.transmit()
    
    def stopTransmit(self):
        self.done = True

    def buttonClickedEvent(self):
        if self.done:
            self.done = False
            self.runTransmit()
            self.buttonWidget.setText("Stop")
            
        elif not self.done:
            sdr = app.sdr
            sdr.tx_destroy_buffer()
            self.stopTransmit()
            self.buttonWidget.setText("Transmit")
            
    def _initFrequencyAddingWidget(self):
        self.tableWidget = QtWidgets.QTableWidget()
        self.buttonWidget = QtWidgets.QPushButton()
        
        self.tableWidget.cellChanged.connect(self.changeFrequencyArray)
        self.tableWidget.setRowCount(3)
        self.tableWidget.setColumnCount(1)
        self.tableWidget.resize(self.tableWidget.sizeHint())
        
        self.gainWidget = QLabeledSlider(Qt.Horizontal)
        self.gainWidget.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.gainWidget.setMinimum(-90)
        self.gainWidget.setMaximum(-10)
        self.gainWidget.setValue(-50)
        self.gainWidget.setSingleStep(1)
        self.gainWidget.valueChanged.connect(self.gainChange)
        
        self.buttonWidget.setText("Transmit")
        self.buttonWidget.clicked.connect(self.buttonClickedEvent)
        
        
        frequencyAddingBox = QtWidgets.QWidget()
        frequencyAddingLayout = QtWidgets.QGridLayout()
        
        frequencyAddingLayout.addWidget(self.tableWidget,0,0)
        frequencyAddingLayout.addWidget(self.buttonWidget,1,0)
        frequencyAddingLayout.addWidget(self.gainWidget,2,0)
        
        
        frequencyAddingBox.setLayout(frequencyAddingLayout)
        
        self.mainLayout.addWidget(frequencyAddingBox)
        
    def changeFrequencyArray(self,row,column):
        self.frequencyTable[row,column] = self.tableWidget.item(row,column).text()
        
        self.signal = np.zeros((self.buffer))
        for row in self.frequencyTable:
            f = row[0]
            w = 2 * np.pi * f
            data = np.sin(w*self.bufferAranged*0.001)
            self.signal += data
        
        self.signal = self.normalize(self.signal)
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.plot(self.bufferAranged,self.signal, '*-')
        self.canvas.draw()
        #self.ax.plot(self.signal, '*-')
        #self.canvas.draw()
    
    def transmit(self):
        if app.isPlutoRunning:
            
            sdr = app.sdr
            sdr.tx_rf_bandwidth = int(app.sampleRate) # filter cutoff, just set it to the same as sample rate
            sdr.tx_lo = int(app.center_freq)
            if self.gain > 0:
                return -1
            sdr.tx_hardwaregain_chan0 = self.gain # Increase to increase tx power, valid range is -90 to 0 dB

            N = 1000 # number of samples to transmit at once
            t = np.arange(N)/app.sampleRate
            samples = None
            for row in self.frequencyTable:
                if samples is None and row[0] is not None:
                    samples = 0.5*np.exp(2.0j*np.pi*row[0]*1e6*t)
                elif row[0] is not None:
                    samples += 0.5*np.exp(2.0j*np.pi*row[0]*1e6*t) # Simulate a sinusoid of 100 kHz, so it should show up at 915.1 MHz at the receiver
            
            samples = self.normalize(samples)
            print(samples)
            samples *= 2**14 # The PlutoSDR expects samples to be between -2^14 and +2^14, not -1 and +1 like some SDRs
            
            sdr.tx_cyclic_buffer = True # Enable cyclic buffers
            print("transmiting!")
            sdr.tx(samples)
            print("transmited!")
            
        else:
            print("Pluto not running!")

    def closeEvent(self, event):
        app.close()
        event.accept()


if __name__ == "__main__":
    # Check whether there is already a running QApplication (e.g., if running
    # from an IDE).
    qapp = QtWidgets.QApplication.instance()
    if not qapp:
        qapp = QtWidgets.QApplication(sys.argv)

    app = ApplicationWindow()
    plot = TransmitWindow()
    app.show()
    plot.show()
    app.activateWindow()
    app.raise_()
    qapp.exec()