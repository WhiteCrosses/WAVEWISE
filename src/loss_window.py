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

    def __init__(self, parent, main_window, qapp):
        super().__init__()
        self.qapp = qapp
        self.parent = parent
        self.result = np.empty(0)
        self.main_window = main_window
        self.currFreq = 0
        self.createWidgets()

    def createWidgets(self):
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

    def transmit(self, sdr, selectedFreq, gain, t):
        print("running")

        base_raw = sdr.rx()

        self.freqs, self.data = signal.periodogram(
            base_raw, self.parent.app.sampleRate)
        self.data = np.where(self.data > 0.00000000001, self.data, -10)
        self.data = 10 * np.log10(np.abs(self.data)**2)

        self.freqs = self.freqs * 1e6 + 1e6

        base = self.find_nearest(self.freqs, freq * 1e6)
        base_start_idx = np.where(self.freqs == base - 1e6*(1/3))[0][0]
        base_stop_idx = np.where(self.freqs == base + 1e6*(1/3))[0][0]

        base_data = self.data[base_start_idx:base_stop_idx]
        base_data = np.max(base_data)

        self.samples = 0.5*np.exp(2.0j*np.pi*freq*1e6*t)
        self.samples = self.parent.normalize(self.samples)
        self.samples *= 2**14

        sdr.tx_cyclic_buffer = True
        sdr.tx(self.samples)
        for x in range(0, 10):
            raw_data = sdr.rx()
        self.signal = sdr.rx()
        sdr.tx_destroy_buffer()

        self.freqs, self.data = signal.periodogram(
            self.signal, self.parent.app.sampleRate)
        self.data = np.where(self.data > 0.00000000001, self.data, -10)
        self.data = 10 * np.log10(np.abs(self.data)**2)

        self.freqs = self.freqs + freq * 1e6 + 1e6

        peaks = signal.find_peaks(self.data, height=-40)

        transmitted = self.find_nearest(self.freqs, freq * 1e6)
        transmitted_idx = np.where(self.freqs == transmitted)
        transmitted_data = self.data[transmitted_idx]

        print(f"base: {base_data},found: {transmitted_data}")
        print(f"loss: {transmitted_data - base_data}")
        print("=====================================")
        return base_data, transmitted_data
        # append to self.result value of freq selected

    def run(self):
        self.currFreq = self.selectedRange[0]

        sdr = self.parent.app.sdr
        N = 1024
        t = np.arange(N)/self.parent.app.sampleRate

        freq_array = np.arange(
            self.selectedRange[0]*1e6, self.selectedRange[1]*1e6, self.stepSelector.value()*1e6)

        self.result = np.empty_like(freq_array)
        idx = 0

        for i in range(self.selectedRange[0], self.selectedRange[1], self.stepSelector.value()):
            base, val = self.transmit(
                sdr, i, self.gainSelector.value(), t)
            val = val - base
            print(val, base)
            self.result[idx] = val
            idx += 1
            self.figure.clear()
            self.ax = self.figure.add_subplot(111)
            self.ax.plot(freq_array, self.result)
            self.canvas.draw()
            # check if canvas completed drawing
            time.sleep(0.1)

    def rangeChangeSlider(self):
        self.startFreqBox.setValue(self.rangeSelectorSlider.value()[0])
        self.endFreqBox.setValue(self.rangeSelectorSlider.value()[1])
        self.selectedRange = self.rangeSelectorSlider.value()

    def rangeChangeBox(self):
        self.selectedRange = (self.startFreqBox.value(),
                              self.endFreqBox.value())
        self.rangeSelectorSlider.setValue(self.selectedRange)
