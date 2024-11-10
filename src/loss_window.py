from superqt import QLabeledRangeSlider, QLabeledSlider  # type: ignore
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvas
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt  # type: ignore
import numpy as np
import threading
import time
from scipy import signal


# TODO detecting sent signal works. it prints value of received signal. now iteration through range of frequencies and displaying them

class LossWindow(QtWidgets.QWidget):
    """!
    @brief [Description de la classe]

    ## Héritage : 
        - Implémente QtWidgets.QWidget => [description]

    """

    def __init__(self, parent, main_window):
        """!
        @brief [Description de la fonction]

        Paramètres : 
            @param self => [description]
            @param parent => [description]
            @param main_window => [description]

        """
        super().__init__()

        self.parent = parent
        self.result = None
        self.main_window = main_window
        self.currFreq = 0
        self.createWidgets()

    def createWidgets(self):
        """!
        @brief [Description de la fonction]

        Paramètres : 
            @param self => [description]

        """
        self.mainLayout = QtWidgets.QFormLayout()

        self.rangeSelectorSlider = QLabeledRangeSlider(Qt.Horizontal)
        self.rangeSelectorSlider.setMinimum(70)
        self.rangeSelectorSlider.setMaximum(6000)
        self.rangeSelectorSlider.setSingleStep(1)
        self.rangeSelectorSlider.valueChanged.connect(self.rangeChangeSlider)

        self.startFreqBox = QtWidgets.QDoubleSpinBox()
        self.startFreqBox.setMinimum(97)
        self.startFreqBox.setMaximum(6000)
        self.startFreqBox.valueChanged.connect(self.rangeChangeBox)

        self.endFreqBox = QtWidgets.QDoubleSpinBox()
        self.endFreqBox.setMinimum(70)
        self.endFreqBox.setMaximum(6000)
        self.endFreqBox.valueChanged.connect(self.rangeChangeBox)

        self.rangeSelectorLayout = QtWidgets.QHBoxLayout()
        self.rangeSelectorLayout.addWidget(self.rangeSelectorSlider)
        self.rangeSelectorLayout.addWidget(self.startFreqBox)
        self.rangeSelectorLayout.addWidget(self.endFreqBox)

        self.rangeSelectorBox = QtWidgets.QWidget()
        self.rangeSelectorBox.setLayout(self.rangeSelectorLayout)

        self.stepSelector = QLabeledSlider(Qt.Horizontal)
        self.stepCounter = QtWidgets.QLabel("0")
        self.gainSelector = QLabeledSlider(Qt.Horizontal)

        self.runButton = QtWidgets.QPushButton("Run!")
        self.runButton.clicked.connect(self.run)

        self.startFreqUnitBox = QtWidgets.QComboBox()
        self.endFreqUnitBox = QtWidgets.QComboBox()

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)

        self.mainLayout.addRow(self.canvas)
        self.mainLayout.addRow("Range to scan [MHz]", self.rangeSelectorBox)
        self.mainLayout.addRow("Single step size [MHz]", self.stepSelector)
        self.mainLayout.addRow("Select gain [mdB]", self.gainSelector)
        self.mainLayout.addRow(self.runButton)
        self.transmitDelay = 1000
        self.setLayout(self.mainLayout)

    def transmitLoop(self):
        # self.runButton.setText("Running...")
        # print("running")
        # self.parent.selectedInex = 0
        # self.parent.gain = self.gainSelector.value()
        # self.parent.frequencyTable[0] = self.selectedRange[0]

        # self.parent.buttonClickedEvent()
        pass

    def transmitPure(self):
        sdr = self.parent.app.sdr
        sdr.tx_cyclic_buffer = True

        threading.Timer(self.transmitDelay, self.stopTransmit).start()

    def stopTransmit(self):
        sdr = self.parent.app.sdr
        sdr.tx_destroy_buffer()

    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def transmit(self):
        print("running")

        sdr = self.parent.app.sdr
        N = 1024
        t = np.arange(N)/self.parent.app.sampleRate
        self.samples = 0.5*np.exp(2.0j*np.pi*self.selectedRange[0]*1e6*t)
        self.samples = self.parent.normalize(self.samples)
        self.samples *= 2**14

        print("transmiting!")
        sdr.tx_cyclic_buffer = True
        sdr.tx(self.samples)
        for x in range(0, 10):
            raw_data = sdr.rx()
        print("transmited!")
        self.signal = sdr.rx()
        sdr.tx_destroy_buffer()

        self.freq, self.data = signal.periodogram(
            self.signal, self.parent.app.sampleRate)
        self.data = np.where(self.data > 0.00000000001, self.data, -10)
        self.data = 10 * np.log10(np.abs(self.data)**2)

        self.freq = self.freq + self.selectedRange[0] * 1e6 + 1e6

        peaks = signal.find_peaks(self.data, height=-40)

        # convert value from 0 to 1000 to range of scan

        print(peaks[0])

        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.plot(self.freq, self.data)

        for i in range(len(peaks[0])):
            self.ax.plot(self.freq[peaks[0][i]], self.data[peaks[0][i]], 'ro')

        self.canvas.draw()

        fndfreq = self.find_nearest(peaks[0], self.selectedRange[0] * 1e6)
        print(self.data[fndfreq])
        # append to self.result value of freq selected

    def run(self):
        self.currFreq = self.selectedRange[0]

        self.transmit()

    def rangeChangeSlider(self):
        self.startFreqBox.setValue(self.rangeSelectorSlider.value()[0])
        self.endFreqBox.setValue(self.rangeSelectorSlider.value()[1])
        self.selectedRange = self.rangeSelectorSlider.value()

    def rangeChangeBox(self):
        self.selectedRange = (self.startFreqBox.value(),
                              self.endFreqBox.value())
        self.rangeSelectorSlider.setValue(self.selectedRange)
