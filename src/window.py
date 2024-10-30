import pyqtgraph as pg  # type: ignore

import adi  # type: ignore

import numpy as np
# Krzywe transmisyjne

from superqt import QLabeledRangeSlider, QDoubleSlider  # type: ignore

from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt, QEvent  # type: ignore

import time
import os
from peak import Peak
import csv

from scipy import signal


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
        self.isSubstracted = False
        self.selectedFreqRange = (70e6, 1000e6)
        self.recType = 0
        self.peakArray = []
        self.delay = 200
        self.selectedRange = "part"
        self.readyData = [[], []]
        self.dataArray = []

        self.press = None
        self.constantPart = 0
        self.sampleRate = int(10e6)
        self.center_freq = int(100e6)
        self.bufferSize = 1024
        self.startFreq = self.center_freq-(self.sampleRate/2)
        self.neededIterations = (
            self.selectedFreqRange[1]-self.selectedFreqRange[0])/self.sampleRate*2 + 1
        self.recMarkerValues = (0, 0)

        self.isClicked = False

        self.desired_level = -70
        self.isPlutoRunning = False
        self.frameTime = 100

        self.isRecording = False
        self.isFirstIteration = False

        self.filterFrame = 5
        self.filterFrameLength = 20

        self.colorMeshMin = self.histMin
        self.colorMeshMax = self.histMax

        self.bgColor = '#0b213b'
        self.setGeometry(0, 0, 1600, 900)

        # self.setFixedWidth(1600)
        # self.setFixedHeight(900)

        self.peaks = []
        self.plot = None
        self._wave_ax_ylim = [-120, 0]
        self._wave_ax_xlim = [70e6, 100e6]

        # self.setStyleSheet("background-color: #0b213b;")

        self.peakFinderInit()
        self.mainWidgetsInit()

        self.setCentralWidget(self._main)

        self.createWidgets()
        self.firstIteration()

    def setPlot(self, plot):
        self.plot = plot

    def firstIteration(self):

        self.x = np.arange(self.startFreq, self.bufferSize+self.startFreq)
        self.signal = self.x * 0
        # change to variables later
        self.imageArray = np.zeros((100, self.bufferSize))
        self.imageArray -= 100

        self._line = self._waterfall_ax.pcolorfast(np.reshape(
            self.imageArray, (100, self.bufferSize)), vmin=self.colorMeshMin, vmax=self.colorMeshMax)
        self._waterfall_ax.invert_yaxis()
        ticks = np.array([0, 1, 2, 3, 4])
        _ticks = ticks/4*self.sampleRate+self.startFreq

        self.calculateWaterfallNodes()
        # self._waterfall_ax.xaxis.set_major_locator(ticker.FixedLocator(_ticks))

        self._wave_ax.set_ylim(self._wave_ax_ylim)
        self.ocid = self._wave_ax.figure.canvas.mpl_connect(
            'button_press_event', self.onClick)
        self.orid = self._wave_ax.figure.canvas.mpl_connect(
            'button_release_event', self.onRelease)
        self.movid = self._wave_ax.figure.canvas.mpl_connect(
            'motion_notify_event', self.onMotion)
        self.scrlid = self._wave_ax.figure.canvas.mpl_connect(
            'scroll_event', self.onScroll)

    def onClick(self, event):
        self.isClicked = True
        self.press = (event.xdata, event.ydata)

    def onRelease(self, event):
        self.isClicked = False

    def onMotion(self, event):
        if self.isClicked:
            xpress, ypress = self.press
            dx = event.xdata - xpress
            dy = event.ydata - ypress

            self._wave_ax_ylim[0] = self._wave_ax_ylim[0] - dy
            self._wave_ax_ylim[1] = self._wave_ax_ylim[1] - dy

            self._wave_ax_xlim[0] = self._wave_ax_xlim[0] - dx
            self._wave_ax_xlim[1] = self._wave_ax_xlim[1] - dx

            self._wave_ax.set_ylim(self._wave_ax_ylim)
            self._wave_ax.set_xlim(self._wave_ax_xlim)

            self._wave_ax.figure.canvas.draw()

    def onScroll(self, event):
        increment = 1 if event.button == 'up' else -1
        # self._wave_ax_ylim[0] = self._wave_ax_ylim[0] - increment * 10
        # self._wave_ax_ylim[1] = self._wave_ax_ylim[1] + increment * 10

        self._wave_ax_xlim[0] = self._wave_ax_xlim[0] - increment * 1e6
        self._wave_ax_xlim[1] = self._wave_ax_xlim[1] + increment * 1e6

        self._wave_ax.set_ylim(self._wave_ax_ylim)
        self._wave_ax.set_xlim(self._wave_ax_xlim)

        self._wave_ax.figure.canvas.draw()

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
        # filter cutoff, just set it to the same as sample rate
        self.sdr.rx_rf_bandwidth = self.sampleRate
        self.sdr.rx_lo = self.center_freq
        self.sdr.gain_control_mode_chan0 = "manual"  # turn off AGC
        gain = 50.0  # allowable range is 0 to 74.5 dB
        self.sdr.rx_hardwaregain_chan0 = gain  # set receive gain
        # this is the buffer the Pluto uses to buffer samples
        self.sdr.rx_buffer_size = self.bufferSize
        self.sdr.rx_rf_bandwidth = int(self.sampleRate)

    def _update_canvas(self):

        if self.isFirstIteration:
            self.plutoInit()
            self.isFirstIteration = False

        if self.isPlutoRunning:
            if self.selectedRange == "part":
                self._wave_ax.clear()
                self._wave_ax_xlim = (
                    self.center_freq-self.sampleRate/2, self.center_freq+self.sampleRate/2)

                self._wave_ax.set_xlim(self._wave_ax_xlim)

                self.getData()

                self._wave_ax.set_ylabel("dBm", color='white')
                self._wave_ax.set_ylim(self._wave_ax_ylim)

                def millions(x, pos):
                    if x < 1e6:
                        return '%1.1fkHz' % (x*1e-3)
                    elif x >= 1e6 and x < 1e9:
                        return '%1.1fMHz' % (x*1e-6)
                    elif x >= 1e9:
                        return '%1.4fGHz' % (x*1e-9)

                formatter = FuncFormatter(millions)
                self._wave_ax.xaxis.set_major_formatter(formatter)

                dataarray = list(zip(self.freq, self.data))
                dataarray = sorted(dataarray, key=lambda x: x[0])

                self.freq, self.data = zip(*dataarray)

                self._wave_ax.plot(self.freq, self.data, color='blue',)

                self.renderPeaks()

                self._line.set_array(np.reshape(
                    self.imageArray, (100, self.bufferSize)))
                self._line.set(clim=(self.colorMeshMin, self.colorMeshMax))

                self._wave_ax.figure.canvas.draw()
                self._line.figure.canvas.draw()

                for row in range(10):
                    self.pRapTable.setItem(
                        row, 0, QtWidgets.QTableWidgetItem(""))
                    self.pRapTable.setItem(
                        row, 1, QtWidgets.QTableWidgetItem(""))
                    self.pRapTable.setItem(
                        row, 2, QtWidgets.QTableWidgetItem(""))

                self.peakArray.sort(key=lambda x: x.frequency)

                if len(self.peakArray) > 1:
                    for x in range(len(self.peakArray)-1):
                        self.peakArray[x+1].distanceToNext = round(
                            (self.freq[self.peakArray[x+1].frequency] - self.freq[self.peakArray[x].frequency])/1e6, 2)

                freqItem = QtWidgets.QTableWidgetItem("Frequency")
                freqItem.setTextAlignment(Qt.AlignHCenter)
                powItem = QtWidgets.QTableWidgetItem("Power")
                powItem.setTextAlignment(Qt.AlignHCenter)
                distItem = QtWidgets.QTableWidgetItem("Distance")
                distItem.setTextAlignment(Qt.AlignHCenter)

                self.pRapTable.setItem(0, 0, freqItem)
                self.pRapTable.setItem(0, 1, powItem)
                self.pRapTable.setItem(0, 2, distItem)

                for row in range(len(self.peakArray)):
                    self.pRapTable.setItem(row+1, 0, QtWidgets.QTableWidgetItem(
                        str(round(self.freq[self.peakArray[row].frequency]/1e6, 2))+" MHz"))
                    self.pRapTable.setItem(
                        row+1, 1, QtWidgets.QTableWidgetItem(str(round(self.peakArray[row].power, 2))+" dBm"))
                    self.pRapTable.setItem(
                        row+1, 2, QtWidgets.QTableWidgetItem(str(self.peakArray[row].distanceToNext)+" MHz"))

            else:
                if self.neededIterations > 0:
                    print(f"{self.neededIterations}")
                    # get data and append it to the array
                    self.getData()

                    dataarray = list(zip(self.freq, self.data))
                    dataarray = sorted(dataarray, key=lambda x: x[0])
                    self.freq, self.data = zip(*dataarray)
                    # self.data = np.subtract(self.data, min(self.data)-self.desired_level)

                    self.readyData = np.concatenate(
                        (self.readyData, np.vstack((self.freq, self.data))), axis=1)
                    # self.readyData.extend([self.freq, self.data])

                    self.center_freq += int(self.sampleRate)
                    self.sdr.rx_lo = self.center_freq
                    self.getFrequencyArray()
                    self.neededIterations -= 1

                else:
                    # data ready to show, show it :)
                    self._wave_ax.clear()
                    # self._wave_ax.set_ylim(self._wave_ax_ylim)
                   # self._wave_ax.set_xlim(self._wave_ax_xlim)

                    self._wave_ax.set_ylabel("dBm", color='white')

                    def millions(x, pos):
                        if x < 1e6:
                            return '%1.1fkHz' % (x*1e-3)
                        elif x >= 1e6 and x < 1e9:
                            return '%1.1fMHz' % (x*1e-6)
                        elif x >= 1e9:
                            return '%1.4fGHz' % (x*1e-9)

                    formatter = FuncFormatter(millions)
                    self._wave_ax.xaxis.set_major_formatter(formatter)
                    self.readyData[self.readyData[:, 0].argsort()]
                    self._wave_ax.plot(
                        self.readyData[0, :], self.readyData[1, :], color='blue')
                    self._wave_ax.set_xlim(
                        (min(self.readyData[0, :]), max(self.readyData[0, :])))

                    # self.renderPeaks()

                    self._wave_ax.figure.canvas.draw()

                    # for row in range(10):
                    # self.pRapTable.setItem(row,0,QtWidgets.QTableWidgetItem(""))

                    # for row in range(len(self.peakArray)):
                    # self.pRapTable.setItem(row,0,QtWidgets.QTableWidgetItem(str(self.freq[self.peakArray[row].frequency])))
                    self.readyData = [[], []]
                    self.center_freq = int(
                        self.selectedFreqRange[0]+self.sampleRate/2)
                    self.neededIterations = (
                        self.selectedFreqRange[1]-self.selectedFreqRange[0])/self.sampleRate*2 + 1

    def renderPeaks(self):

        x, y = self.checkPeaks()

        # Iterate through peaks to assign them new positions
        m = 0
        for i in x:
            found = False
            for j in self.peakArray:
                if j.min <= i <= j.max:
                    j.found(i, y['peak_heights'][m])
                    found = True
                    break

            if not found:
                self.peakArray.append(
                    Peak(i, y['peak_heights'][m], self.pSetDistance))
            m += 1

        for j in self.peakArray:
            if not j.isChecked:
                j.tickCounter += 1
                if j.tickCounter > 10:
                    self.peakArray.remove(j)

        for j in self.peakArray:
            j.isChecked = False

        for i in self.peakArray:
            self._wave_ax.axvline(self.freq[i.frequency], color=i.color)

    def getFrequencyArray(self):
        self.freq = np.fft.fftfreq(1024, d=1/self.sampleRate)
        self.freq = self.freq + self.center_freq

    def getData(self):
        if self.isPlutoRunning:
            self.signal = self.sdr.rx()
            self.processData2()
            # self.getFrequencyArray()

    def dataFilter(self):
        kernel = np.ones(10)/10
        self.data = np.convolve(self.data, kernel, mode='same')

    def processData2(self):

        self.freq, self.data = signal.periodogram(self.signal, self.sampleRate)

        self.freq = self.freq + self.center_freq

        self.data = np.where(self.data > 0.00000000001, self.data, -10)
        self.data = 10 * np.log10(np.abs(self.data)**2)

        thereshold = np.average(self.data)*0.9
        data_tmp = self.data
        for i in range(len(self.data)):
            if self.data[i] < thereshold:
                data_tmp[i] = np.average(self.data)  # * 0.9
            else:
                data_tmp[i] = self.data[i]

        self.data = data_tmp
        # Smooth out noise
        self.dataFilter()

        self.dataFilter()
        self.addToImageArray()

    def processData(self):

        # process pulled data from pluto
        # Fast Fourier transform
        # self.data = fft(self.signal)
        # self.signal = self.signal - np.mean(self.signal)

        # self.data = np.abs(np.fft.fft(self.signal))**2 / (self.bufferSize*self.sampleRate)
        window = np.hamming(self.bufferSize)
        self.signal = self.signal * window
        ft = np.fft.fft(self.signal)
        self.data = np.roll(ft, int(self.bufferSize//2))

        # print(type(self.data))
        # self.data = (ft.real**2 + ft.imag**2)/ (self.bufferSize*self.sampleRate)

        # Converting to logarythmic scale
        # self.data = (20 * np.log10(self.data))
        # self.data = 2/self.sampleRate*self.bufferSize*np.abs(self.data)
        # print(self.data)
        # average out the sudden drops in power
        # If value is below 0.9 average it is set to average
        thereshold = np.average(self.data)*0.9
        data_tmp = self.data
        for i in range(len(self.data)):
            if self.data[i] < thereshold:
                data_tmp[i] = np.average(self.data)  # * 0.9
            else:
                data_tmp[i] = self.data[i]

        self.data = data_tmp
        # Smooth out noise
        self.dataFilter()

        self.addToImageArray()

    def addToImageArray(self):
        self.imageArray = np.roll(self.imageArray, -1, axis=0)
        self.imageArray[0, :] = np.roll(
            self.data, int(self.data.size/2), axis=0)

    def checkPeaks(self):
        x, y = signal.find_peaks(
            self.data,
            height=self.pSetHeight,
            threshold=self.pSetThereshold,
            distance=self.pSetDistance,
            # prominence=self.pSetProminence,
            width=self.pSetWidth,
            # wlen=self.pSetWlen,
            rel_height=self.pSetRelHeight,
            plateau_size=self.pSetPlateauSize)
        return x, y

    def calculateWaterfallNodes(self):

        def millions(x):
            if x < 1e6:
                return '%1.1fkHz' % (x*1e-3)
            elif x >= 1e6 and x < 1e9:
                return '%1.1fMHz' % (x*1e-6)
            elif x >= 1e9:
                return '%1.4fGHz' % (x*1e-9)

        ticks = np.array([0, 1, 2, 3, 4])
        _ticks = ticks/4*self.sampleRate+self.startFreq

        formattedTicks = ['', '', '', '', '']
        for i in range(len(_ticks)):
            formattedTicks[i] = millions(_ticks[i])

        tickLoc = np.arange(0, 1024, 1024/4)
        tickLoc = np.append(tickLoc, 1024)
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
        if value <= 0:
            pass
        else:
            self.center_freq = int(value * 1e6)
            self.startFreq = self.center_freq-(self.sampleRate/2)
            self.sdr.rx_lo = self.center_freq

        self.calculateWaterfallNodes()

        self.getFrequencyArray()

    def plutoStartChangeState(self):
        if self.isPlutoRunning:
            self.isPlutoRunning = False

            self.mainWidgetsInit()
            self.setCentralWidget(self._main)
            self.createWidgets()
            self.firstIteration()
        elif not self.isPlutoRunning:
            self.isPlutoRunning = True
            self.plutoInit()
            self.firstIteration()
            self.plutoStartButton.setText("Stop Pluto!")

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
        self.pSetDistance = int(
            (self.bufferSize/self.sampleRate)*1e6*self.pSetDistanceSlider.value())
        if self.pSetDistance < 1:
            self.pSetDistance = 1
        self.pSetDistanceSpinBox.setValue(self.pSetDistanceSlider.value())

    def pSetDistanceChangeSpinBox(self):
        self.pSetDistance = int(
            (self.bufferSize/self.sampleRate)*1e6*self.pSetDistanceSpinBox.value())
        if self.pSetDistance < 1:
            self.pSetDistance = 1
        self.pSetDistanceSlider.setValue(self.pSetDistanceSpinBox.value())

    def pSetTheresholdChange(self):
        self.pSetThereshold = -self.bufferSize / \
            self.sampleRate*self.pSetTheresholdSlider.value()
        self.pSetTheresholdSpinBox.setValue(self.pSetTheresholdSlider.value())

    def pSetTheresholdChangeSpinBox(self):
        self.pSetThereshold = -self.bufferSize / \
            self.sampleRate*self.pSetTheresholdSpinBox.value()
        self.pSetTheresholdSlider.setValue(self.pSetTheresholdSpinBox.value())

    def pSetProminenceChange(self):
        self.pSetProminence = self.pSetProminenceSlider.value()

    def pSetWidthChange(self):
        self.pSetWidth = int(self.bufferSize/self.sampleRate *
                             1e6*self.pSetWidthSpinBox.value())
        self.pSetWidthSpinBox.setValue(self.pSetWidthSlider.value())

    def pSetWidthChangeSpinBox(self):
        self.pSetWidth = int(self.bufferSize/self.sampleRate *
                             1e6*self.pSetWidthSpinBox.value())
        self.pSetWidthSlider.setValue(self.pSetWidthSpinBox.value())

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
            pass

        return False

    def recTypeChange(self, index):
        self.recType = index

    def recRangeChange(self):
        self.recMarkerValues = self.recRangeSelector.value()

    def savePeaks(self, filenameP):
        with open(filenameP, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Frequency [Hz]", "Power [dBm]"])
            m = 0
            for i in self.peakArray:
                writer.writerow([self.freq[i.frequency], i.power])
            file.close()

    def saveWaterfall(self, filenameW):
        self._waterfall_ax.figure.savefig(filenameW)

    def saveData(self, filenameD):
        with open(filenameD, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Frequency [Hz]", "Power [dBm]"])
            m = 0
            for i in range(len(self.data)):
                writer.writerow([self.freq[i], self.data[i]])
            file.close()

    def recResolve(self):

        destFolder = "./output"

        if not os.path.exists(destFolder):
            os.makedirs(destFolder)

        date = time.strftime("%Y_%m_%d_%H_%M_%S")
        filenameP = destFolder+"/"+"PEAKS_"+date+".csv"
        filenameW = destFolder+"/"+date+".png"
        filenameD = destFolder+"/"+date+".csv"

        match self.recType:
            case 0:
                self.savePeaks(filenameP)
            case 1:
                self.saveWaterfall(filenameW)
            case 2:
                self.saveData(filenameD)
            case 3:
                self.savePeaks(filenameP)
                self.saveData(filenameD)
                self.saveWaterfall(filenameW)

    def pRapSelectPeak(self, row, column):
        value = self.pRapTable.item(row, column).text()
        value = float(value)
        for peak in self.peakArray:
            if peak.min <= float(value) and peak.max >= float(value):
                peak.color = 'blue'
                return 0

    def rangeSwitchChange(self):
        if self.rangeSwitchButton.text() == "Part of spectrum":
            self.rangeSwitchButton.setText("Full spectrum")
            self.selectedRange = "full"
        else:
            self.rangeSwitchButton.setText("Part of spectrum")
            self.selectedRange = "part"

        self.mainWidgetsInit()
        self.setCentralWidget(self._main)
        self.createWidgets()
        self.firstIteration()

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

        # Peak settings
        # height
        self.pSetHeightSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetHeightSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetHeightSlider.setTickInterval(5)
        self.pSetHeightSlider.setMinimum(-100)
        self.pSetHeightSlider.setMaximum(0)
        self.pSetHeightSlider.setSingleStep(5)
        self.pSetHeightSlider.valueChanged.connect(self.pSetHeightChange)
        self.pSetHeightSlider.description = "Declare minimum height to peak be detected. Values ranging from -100 to 0."
        self.pSetHeightSlider.installEventFilter(self)

        self.pSetHeightSpinBox = QtWidgets.QSpinBox()
        self.pSetHeightSpinBox.valueChanged.connect(
            self.pSetHeightChangeSpinBox)
        self.pSetHeightSpinBox.setRange(-100, 0)
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

        # thereshold
        self.pSetTheresholdSlider = QDoubleSlider(Qt.Horizontal)
        self.pSetTheresholdSlider.setTickPosition(
            QtWidgets.QSlider.TicksBothSides)
        self.pSetTheresholdSlider.setTickInterval(1.5)
        self.pSetTheresholdSlider.setMinimum(0)
        self.pSetTheresholdSlider.setMaximum(15)
        self.pSetTheresholdSlider.setSingleStep(0.1)
        self.pSetTheresholdSlider.valueChanged.connect(
            self.pSetTheresholdChange)
        self.pSetTheresholdSlider.description = "Declare minimum diference in the values between samples"
        self.pSetTheresholdSlider.installEventFilter(self)

        self.pSetTheresholdSpinBox = QtWidgets.QDoubleSpinBox()
        self.pSetTheresholdSpinBox.valueChanged.connect(
            self.pSetTheresholdChangeSpinBox)
        self.pSetTheresholdSpinBox.setRange(0, 15)
        self.pSetTheresholdSpinBox.setSingleStep(0.1)

        self.pSetTheresholdSlider.setValue(2)
        self.pSetTheresholdSpinBox.setValue(2)

        self.pSetTheresholdBoxLayout = QtWidgets.QHBoxLayout()
        self.pSetTheresholdBoxLayout.addWidget(self.pSetTheresholdSlider)
        self.pSetTheresholdBoxLayout.addWidget(self.pSetTheresholdSpinBox)
        self.pSetTheresholdBox = QtWidgets.QHBoxLayout()
        self.pSetTheresholdBox.addLayout(self.pSetTheresholdBoxLayout)

        # distance
        self.pSetDistanceSlider = QDoubleSlider(Qt.Horizontal)
        # self.pSetDistanceSlider = QtWidgets.QSlider(Qt.Horizontal)
        self.pSetDistanceSlider.setTickPosition(
            QtWidgets.QSlider.TicksBothSides)
        self.pSetDistanceSlider.setTickInterval(1)
        self.pSetDistanceSlider.setMinimum(0)
        self.pSetDistanceSlider.setMaximum(10)
        self.pSetDistanceSlider.setSingleStep(0.1)
        self.pSetDistanceSlider.valueChanged.connect(self.pSetDistanceChange)
        # self.pSetDistanceSlider.setToolTip("Required minimal horizontal distance (>= 1) in samples between neighbouring peaks.\nSmaller peaks are removed first until the condition is fulfilled for all remaining peaks.")
        self.pSetDistanceSlider.description = "Required minimal horizontal distance (>= 1) in samples between neighbouring peaks.\nSmaller peaks are removed first until the condition is fulfilled for all remaining peaks."
        self.pSetDistanceSlider.installEventFilter(self)

        self.pSetDistanceSpinBox = QtWidgets.QDoubleSpinBox()
        self.pSetDistanceSpinBox.valueChanged.connect(
            self.pSetDistanceChangeSpinBox)
        self.pSetDistanceSpinBox.setRange(0, 10)
        self.pSetDistanceSpinBox.setSingleStep(0.1)

        self.pSetDistanceBoxLayout = QtWidgets.QHBoxLayout()
        self.pSetDistanceBoxLayout.addWidget(self.pSetDistanceSlider)
        self.pSetDistanceBoxLayout.addWidget(self.pSetDistanceSpinBox)
        self.pSetDistanceBox = QtWidgets.QHBoxLayout()
        self.pSetDistanceBox.addLayout(self.pSetDistanceBoxLayout)

        # width
        self.pSetWidthSlider = QDoubleSlider(Qt.Horizontal)
        self.pSetWidthSlider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.pSetWidthSlider.setTickInterval(0.1)
        self.pSetWidthSlider.setMinimum(0)
        self.pSetWidthSlider.setMaximum(1)
        self.pSetWidthSlider.setSingleStep(0.01)
        self.pSetWidthSlider.valueChanged.connect(self.pSetWidthChange)
        # self.pSetWidthSlider.setToolTip(f"Required width of peaks in percentage of all samples ({self.bufferSize}).")
        self.pSetWidthSlider.description = f"Required width of peaks in MHz)."
        self.pSetWidthSlider.installEventFilter(self)

        self.pSetWidthSpinBox = QtWidgets.QDoubleSpinBox()
        self.pSetWidthSpinBox.valueChanged.connect(self.pSetWidthChangeSpinBox)
        self.pSetWidthSpinBox.setRange(0, 1)
        self.pSetWidthSpinBox.setSingleStep(0.01)

        self.pSetWidthBoxLayout = QtWidgets.QHBoxLayout()
        self.pSetWidthBoxLayout.addWidget(self.pSetWidthSlider)
        self.pSetWidthBoxLayout.addWidget(self.pSetWidthSpinBox)
        self.pSetWidthBox = QtWidgets.QHBoxLayout()
        self.pSetWidthBox.addLayout(self.pSetWidthBoxLayout)

        self.statusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.statusBar)
        self.statusBar.showMessage(" ")

        self.recTypeChoice = QtWidgets.QComboBox()
        self.recTypeChoice.addItems(("Save peaks",
                                    "Save waterfall",
                                     "Save data frame",
                                     "Save all"))
        self.recTypeChoice.currentIndexChanged.connect(self.recTypeChange)
        self.recTypeChoice.description = "Declare requested data to be saved"
        self.recTypeChoice.installEventFilter(self)

        self.recRangeSelector = QLabeledRangeSlider(Qt.Horizontal)
        self.recRangeSelector.setValue((0, 121))
        self.recRangeSelector.valueChanged.connect(self.recRangeChange)

        # Peak tracking table
        self.pRapTable = QtWidgets.QTableWidget(10, 3)
        self.pRapTable.description = "List of peaks detected."
        self.pRapTable.installEventFilter(self)
        self.pRapTable.cellPressed.connect(self.pRapSelectPeak)

        # Start/Stop button for recording to csv
        self.recButton = QtWidgets.QPushButton("Save!")
        self.recButton.clicked.connect(self.recChangeRecordingState)
        recButtonLayout = QtWidgets.QGridLayout()
        recButtonLayout.addWidget(self.recButton, 0, 0)
        recButtonLayout.addWidget(self.recTypeChoice, 0, 1)
        self.recButton.description = "Save data to the file."
        self.recButton.installEventFilter(self)
        #   recButtonLayout.addWidget(self.recRangeSelector,1,0,1,2)

        recordingBoxWidget = QtWidgets.QGroupBox("Recording")
        recordingBoxWidget.setLayout(recButtonLayout)

        if self.selectedRange == "part":
            self.rangeSwitchButton = QtWidgets.QPushButton("Part of spectrum")
        else:
            self.rangeSwitchButton = QtWidgets.QPushButton("Full spectrum")
        self.rangeSwitchButton.clicked.connect(self.rangeSwitchChange)

        self.plutoStartButton = QtWidgets.QPushButton("Start Pluto!")
        self.plutoStartButton.installEventFilter(self)
        self.plutoStartButton.description = "Start the device"
        self.plutoStartButton.clicked.connect(self.plutoStartChangeState)
        plutoStartButtonLayout = QtWidgets.QHBoxLayout()
        plutoStartButtonLayout.addWidget(self.plutoStartButton)
        plutoStartBoxWidget = QtWidgets.QGroupBox("I/O")
        plutoStartBoxWidget.setLayout(plutoStartButtonLayout)

        layoutpSetBox = QtWidgets.QFormLayout()
        layoutpSetBox.addRow("Height [dBm]", self.pSetHeightBox)
        layoutpSetBox.addRow("Thereshold [dBm]", self.pSetTheresholdBox)
        layoutpSetBox.addRow("Distance [MHz]", self.pSetDistanceBox)
        # layoutpSetBox.addRow("Prominence", self.pSetProminenceSlider)
        layoutpSetBox.addRow("Width [MHz]", self.pSetWidthBox)
        # layoutpSetBox.addRow("Wlen", self.pSetWlenSlider)
        # layoutpSetBox.addRow("Relative Height", self.pSetRelHeightSlider)
        # layoutpSetBox.addRow("Plateau Size", self.pSetPlateauSizeSlider)

        layoutpRapBox = QtWidgets.QTableWidget

        pSetBoxWidget = QtWidgets.QGroupBox("Peak settings")
        pSetBoxWidget.setLayout(layoutpSetBox)

        layoutMainBox = QtWidgets.QHBoxLayout(self._main)
        layoutMiddleBox = QtWidgets.QVBoxLayout()
        layoutLeftBox = QtWidgets.QVBoxLayout()

        wave_canvas = FigureCanvas(Figure(figsize=(10, 6)))
        waterfall_canvas = FigureCanvas(Figure(figsize=(10, 6)))

        layoutMiddleBox.addWidget(wave_canvas)
        if self.selectedRange == "part":
            layoutMiddleBox.addWidget(waterfall_canvas)

        layoutSamplingSetBox = QtWidgets.QFormLayout()
        layoutSamplingSetBox.addRow("Sample rate [MHz]", self.sampleRateBox)
        layoutSamplingSetBox.addRow(
            "Center Frequency [MHz]", self.centerFreqBox)
        samplingSetBoxWidget = QtWidgets.QGroupBox("Peak settings")
        samplingSetBoxWidget.setLayout(layoutSamplingSetBox)

        layoutLeftBox.addWidget(plutoStartBoxWidget)
        layoutLeftBox.addWidget(recordingBoxWidget)
        # layoutLeftBox.addWidget(self.rangeSwitchButton)
        if self.selectedRange == "part":
            layoutLeftBox.addWidget(samplingSetBoxWidget)
            layoutLeftBox.addWidget(pSetBoxWidget)
            layoutLeftBox.addWidget(self.pRapTable)

        self._leftBox.setLayout(layoutLeftBox)
        self._middleBox.setLayout(layoutMiddleBox)

        layoutMainBox.addWidget(self._leftBox)
        layoutMainBox.addWidget(self._middleBox)
        if self.selectedRange == "part":
            layoutMainBox.addWidget(self.hist)

        # button1.clicked.connect(self.moveCenter)

        wave_canvas.figure.set_facecolor(self.bgColor)
        waterfall_canvas.figure.set_facecolor(self.bgColor)

        with plt.rc_context({'axes.edgecolor': 'white', 'xtick.color': 'white', 'ytick.color': 'white', 'figure.facecolor': self.bgColor}):

            self._wave_ax = wave_canvas.figure.subplots()

            self._waterfall_ax = waterfall_canvas.figure.subplots()
            self._wave_ax.tick_params(colors='white', which='both')

        self._timer = waterfall_canvas.new_timer(self.delay)
        self._timer.add_callback(self._update_canvas)
        self._timer.start()

    def closeEvent(self, event):
        self.plot.close()
        event.accept()
