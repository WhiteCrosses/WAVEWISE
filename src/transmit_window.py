import numpy as np
# Krzywe transmisyjne

from superqt import QLabeledSlider  # type: ignore

from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvas
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt  # type: ignore
import threading


class TransmitWindow(QtWidgets.QWidget):
    def __init__(self, app):
        super().__init__()
        self.app = app
        self.frequencyTable = np.array(np.zeros((3, 2)))
        self.buffer = 1024
        self.gain = -50
        self.selectedIndex = 2
        self.bufferAranged = np.array(np.arange(0, 1024, 1))
        self.mainLayout = QtWidgets.QVBoxLayout()
        self.done = True
        self.cyclicBreaker = False
        self.transmitDelay = 2
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

    def normalize(self, arr):
        # Find the maximum absolute magnitude in the array
        max_magnitude = np.max(np.abs(arr))
        # Divide the array by the maximum magnitude
        normalized_arr = arr / max_magnitude
        return normalized_arr

    def buttonClickedEvent(self):
        if self.buttonWidget.text() == "Transmit":
            self.transmit()
            if self.selectedIndex == 0:
                sdr = self.app.sdr
                sdr.tx_destroy_buffer()
            else:
                self.buttonWidget.setText("Stop")

        else:
            if self.selectedIndex == 1:
                self.cyclicBreaker = True
                # self.transmitCallback()
            sdr = self.app.sdr
            sdr.tx_destroy_buffer()
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

        self.cyclicSelector = QtWidgets.QComboBox()
        self.cyclicSelector.addItems(("Pulse",
                                      "Cyclic",
                                      "Continuous"))

        self.cyclicSelector.currentIndexChanged.connect(self.indexChanged)

        self.cyclicSelectorLabel = QtWidgets.QLabel()
        self.cyclicSelectorLabel.setText("Select transmission type:")

        frequencyAddingBox = QtWidgets.QWidget()
        frequencyAddingLayout = QtWidgets.QGridLayout()

        frequencyAddingLayout.addWidget(self.tableWidget, 0, 0, 1, 2)
        frequencyAddingLayout.addWidget(self.cyclicSelectorLabel, 1, 0)
        frequencyAddingLayout.addWidget(self.cyclicSelector, 1, 1)
        frequencyAddingLayout.addWidget(self.buttonWidget, 2, 0, 1, 2)
        frequencyAddingLayout.addWidget(self.gainWidget, 3, 0, 1, 2)

        frequencyAddingBox.setLayout(frequencyAddingLayout)

        self.mainLayout.addWidget(frequencyAddingBox)

    def indexChanged(self, index):
        print(f"index changed to: {index}")
        self.selectedIndex = index

    def changeFrequencyArray(self, row, column):
        self.frequencyTable[row, column] = self.tableWidget.item(
            row, column).text()

        self.signal = np.zeros((self.buffer))
        for row in self.frequencyTable:
            f = row[0]
            w = 2 * np.pi * f
            data = np.sin(w*self.bufferAranged*0.001)
            self.signal += data

        # self.signal = self.normalize(self.signal)
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.plot(self.bufferAranged, self.signal, '*-')
        self.canvas.draw()
        # self.ax.plot(self.signal, '*-')
        # self.canvas.draw()

    def transmitCallback(self):
        if self.selectedIndex == 1 and not self.cyclicBreaker:
            print("Adding next recursion")
            self._timer = threading.Timer(
                self.transmitDelay, self.transmitCallback).start()

        print("entered callback")
        sdr = self.app.sdr

        # filter cutoff, just set it to the same as sample rate
        sdr.tx_rf_bandwidth = int(self.app.sampleRate)
        sdr.tx_lo = int(self.app.center_freq)
        if self.gain > 0:
            return -1
        # Increase to increase tx power, valid range is -90 to 0 dB
        sdr.tx_hardwaregain_chan0 = self.gain
        N = 1000  # number of samples to transmit at once
        t = np.arange(N)/self.app.sampleRate
        samples = None
        for row in self.frequencyTable:
            if samples is None and row[0] is not None:
                samples = 0.5*np.exp(2.0j*np.pi*row[0]*1e6*t)
            elif row[0] is not None:
                # Simulate a sinusoid of 100 kHz, so it should show up at 915.1 MHz at the receiver
                samples += 0.5*np.exp(2.0j*np.pi*row[0]*1e6*t)
        samples = self.normalize(samples)

        samples *= 2**14  # The PlutoSDR expects samples to be between -2^14 and +2^14, not -1 and +1 like some SDRs

        sdr.tx_cyclic_buffer = True
        print("transmiting!")
        sdr.tx(samples)
        print("transmited!")
        for x in range(0, 10):
            raw_data = sdr.rx()
        sdr.tx_destroy_buffer()

    def transmit(self):
        if self.app.isPlutoRunning:

            if self.selectedIndex == 0:   # Pulse
                print("Performing pulse transmission")
                self.transmitCallback()

            elif self.selectedIndex == 1:  # Cyclic
                print("Starting cyclic transmission")
                self.cyclicBreaker = False
                self._timer = threading.Timer(
                    self.transmitDelay, self.transmitCallback)
                self._timer.start()

            elif self.selectedIndex == 2:  # Continuous
                print("Starting continuous transmission")
                sdr = self.app.sdr
                # filter cutoff, just set it to the same as sample rate
                sdr.tx_rf_bandwidth = int(self.app.sampleRate)
                sdr.tx_lo = int(self.app.center_freq)
                if self.gain > 0:
                    return -1
                # Increase to increase tx power, valid range is -90 to 0 dB
                sdr.tx_hardwaregain_chan0 = self.gain
                N = 1000  # number of samples to transmit at once
                t = np.arange(N)/self.app.sampleRate
                samples = None
                for row in self.frequencyTable:
                    if samples is None and row[0] is not None:
                        samples = 0.5*np.exp(2.0j*np.pi*row[0]*1e6*t)
                    elif row[0] is not None:
                        # Simulate a sinusoid of 100 kHz, so it should show up at 915.1 MHz at the receiver
                        samples += 0.5*np.exp(2.0j*np.pi*row[0]*1e6*t)
                samples = self.normalize(samples)

                samples *= 2**14  # The PlutoSDR expects samples to be between -2^14 and +2^14, not -1 and +1 like some SDRs

                print("transmiting!")
                sdr.tx_cyclic_buffer = True
                sdr.tx(samples)
                print("transmited!")

        else:
            print("Pluto not running!")

    def closeEvent(self, event):
        self.app.close()
        event.accept()
