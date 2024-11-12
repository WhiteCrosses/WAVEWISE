from superqt import QLabeledRangeSlider, QLabeledSlider  # type: ignore
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvas
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt  # type: ignore
import numpy as np
import threading
import time
from scipy import signal
import adi  # type: ignore

# TODO detecting sent signal works. it prints value of received signal. now iteration through range of frequencies and displaying them


class LossWindow(QtWidgets.QWidget):
    """!
    @brief [Description de la classe]

    ## Héritage :
        - Implémente QtWidgets.QWidget => [description]

    """
    """!
    @brief [Description de la classe]

    ## Héritage :
        - Implémente QtWidgets.QWidget => [description]

    """

    def __init__(self, app):
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]
            @param app => [description]

        """
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]
            @param qapp => [description]

        """
        super().__init__()
        self.app = app
        self.selectedRange = (97, 200)
        self.sampleRate = int(10e6)
        self.center_freq = int(96*1e6)
        self.bufferSize = 2048

        self.plutoInit()
        self.result = np.array([])
        self.currFreq = 0
        self.createWidgets()

    def plutoInit(self):
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]

        """
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]

        """
        self.sdr = adi.Pluto("ip:192.168.2.1")
        self.sdr.sample_rate = self.sampleRate
        # filter cutoff, just set it to the same as sample rate
        self.sdr.rx_rf_bandwidth = self.sampleRate
        self.sdr.rx_lo = self.center_freq
        self.sdr.gain_control_mode_chan0 = "manual"  # turn off AGC
        gain = 50.0  # allowable range is 0 to 74.5 dB
        self.sdr.rx_hardwaregain_chan0 = gain  # set receive gain
        self.sdr.tx_hardwaregain_chan0 = -20
        # this is the buffer the Pluto uses to buffer samples
        self.sdr.rx_buffer_size = self.bufferSize
        self.sdr.rx_rf_bandwidth = int(self.sampleRate)
        self.sampleRate = self.sdr.sample_rate

    def createWidgets(self):
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]

        """
        self.mainLayout = QtWidgets.QFormLayout()

        self.rangeSelectorSlider = QLabeledRangeSlider(Qt.Horizontal)
        self.rangeSelectorSlider.setMinimum(97)
        self.rangeSelectorSlider.setMaximum(6000)
        self.rangeSelectorSlider.setSingleStep(1)
        self.rangeSelectorSlider.setValue((97, 200))
        self.rangeSelectorSlider.valueChanged.connect(self.rangeChangeSlider)

        self.startFreqBox = QtWidgets.QDoubleSpinBox()
        self.startFreqBox.setMinimum(200)
        self.startFreqBox.setMaximum(6000)
        self.startFreqBox.setValue(200)
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
        self.stepSelector.setValue(1)
        self.gainSelector = QLabeledSlider(Qt.Horizontal)
        self.gainSelector.setValue(50)

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

    def find_nearest(self, array, value):
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]
            @param array => [description]
            @param value => [description]

        """
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def normalize(self, arr):
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]
            @param arr => [description]

        """
        # Find the maximum absolute magnitude in the array
        max_magnitude = np.max(np.abs(arr))
        # Divide the array by the maximum magnitude
        normalized_arr = arr / max_magnitude
        return normalized_arr

    def transmit_test(self, sdr, selectedFreq, gain, t, delay):
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]
            @param sdr => [description]
            @param selectedFreq => [description]
            @param gain => [description]
            @param t => [description]

        """
        print("running")

        selectedFreq = int(selectedFreq * 1e6)

        self.samples = 0.5*np.exp(2.0j*np.pi*selectedFreq*t)
        self.samples = self.normalize(self.samples)
        self.samples *= 2**14

        self.sdr.rx_lo = int(selectedFreq - 2e6)

        for x in range(0, 5):
            raw_data = self.sdr.rx()
        base_raw = self.sdr.rx()
        self.sdr.tx_cyclic_buffer = True
        self.sdr.tx(self.samples)

        for x in range(0, delay):
            raw_data = self.sdr.rx()
        sig_raw = self.sdr.rx()
        self.sdr.tx_destroy_buffer()

        self.freqs_base, self.data_base = signal.periodogram(
            base_raw, self.sampleRate)
        self.data_base = np.where(
            self.data_base > 0.00000000001, self.data_base, -10)
        self.data_base = 10 * np.log10(np.abs(self.data_base)**2)

        self.freqs_sig, self.data_sig = signal.periodogram(
            sig_raw, self.sampleRate)
        self.data_sig = np.where(
            self.data_sig > 0.00000000001, self.data_sig, -10)
        self.data_sig = 10 * np.log10(np.abs(self.data_sig)**2)

        self.freqs_sig = self.freqs_sig + selectedFreq - 2e6
        self.freqs_base = self.freqs_base + selectedFreq - 2e6

        base = self.find_nearest(self.freqs_base, selectedFreq)

        base_idx = np.where(self.freqs_base == base)[0][0]
        base_idx += 1
        sigOff = self.data_base[base_idx]
        sigOn = self.data_sig[base_idx]

        print(f"selectedFreq: {selectedFreq}")
        print(f"f: {self.freqs_base[base_idx]}")
        print(f"off: {sigOff},on: {sigOn}")
        print(f"loss: {sigOn - sigOff}")
        print("=====================================")

        return sigOff, sigOn, self.data_base, self.data_sig, self.freqs_base

    def transmit(self, sdr, selectedFreq, gain, t):
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]
            @param sdr => [description]
            @param selectedFreq => [description]
            @param gain => [description]
            @param t => [description]

        """
        print("running")

        selectedFreq = int(selectedFreq * 1e6)

        self.samples = 0.5*np.exp(2.0j*np.pi*selectedFreq*t)
        self.samples = self.normalize(self.samples)
        self.samples *= 2**14

        self.sdr.rx_lo = int(selectedFreq - 2e6)

        for x in range(0, 5):
            raw_data = self.sdr.rx()
        base_raw = self.sdr.rx()
        self.sdr.tx_cyclic_buffer = True
        self.sdr.tx(self.samples)

        for x in range(0, 5):
            raw_data = self.sdr.rx()
        sig_raw = self.sdr.rx()
        self.sdr.tx_destroy_buffer()

        self.freqs_base, self.data_base = signal.periodogram(
            base_raw, self.sampleRate)
        self.data_base = np.where(
            self.data_base > 0.00000000001, self.data_base, -10)
        self.data_base = 10 * np.log10(np.abs(self.data_base)**2)

        self.freqs_sig, self.data_sig = signal.periodogram(
            sig_raw, self.sampleRate)
        self.data_sig = np.where(
            self.data_sig > 0.00000000001, self.data_sig, -10)
        self.data_sig = 10 * np.log10(np.abs(self.data_sig)**2)

        self.freqs_sig = self.freqs_sig + selectedFreq - 2e6
        self.freqs_base = self.freqs_base + selectedFreq - 2e6

        base = self.find_nearest(self.freqs_base, selectedFreq)

        base_idx = np.where(self.freqs_base == base)[0][0]
        base_idx += 1
        sigOff = self.data_base[base_idx]
        sigOn = self.data_sig[base_idx]

        print(f"selectedFreq: {selectedFreq}")
        print(f"f: {self.freqs_base[base_idx]}")
        print(f"off: {sigOff},on: {sigOn}")
        print(f"loss: {sigOn - sigOff}")
        print("=====================================")
        return sigOff, sigOn, self.data_base, self.data_sig, self.freqs_base
        # append to self.result value of freq selected

    def run(self):
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]

        """
        self.currFreq = self.selectedRange[0]
        N = 10000
        t = np.arange(N)/self.sampleRate
        window = np.ones(10)/float(10)
        freq_array = np.arange(
            self.selectedRange[0]*1e6, self.selectedRange[1]*1e6, self.stepSelector.value()*1e6)

        idx = 0
        self.result = np.array([])

        for i in range(self.selectedRange[0], self.selectedRange[1], self.stepSelector.value()):
            sigOff, sigOn, fft_off, fft_on, freqs = self.transmit(
                self.sdr, i, self.gainSelector.value(), t)
            val = sigOn - sigOff - self.gainSelector.value()
            self.result = np.append(self.result, val)
            displayed = np.convolve(self.result, window, 'same')

            idx += 1
            self.figure.clear()
            plt.ion()
            self.ax = self.figure.add_subplot(211)
            self.ax2 = self.figure.add_subplot(212)

            self.ax2.set_xlim(i*1e6 - self.sampleRate/2 - 2e6,
                              i*1e6 + self.sampleRate/2 - 2e6)

            if (idx > 10):
                self.ax.plot(freq_array[0:idx], displayed)

                self.ax2.plot(freqs, fft_on)
                self.ax2.plot(freqs, fft_off)
            self.canvas.draw()
            for x in range(0, 10):
                time.sleep(0.01)

    def rangeChangeSlider(self):
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]

        """
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]

        """
        self.startFreqBox.setValue(self.rangeSelectorSlider.value()[0])
        self.endFreqBox.setValue(self.rangeSelectorSlider.value()[1])
        self.selectedRange = self.rangeSelectorSlider.value()

    def rangeChangeBox(self):
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]

        """
        """!
        @brief [Description de la fonction]

        Paramètres :
            @param self => [description]

        """
        self.selectedRange = (self.startFreqBox.value(),
                              self.endFreqBox.value())
        self.rangeSelectorSlider.setValue(self.selectedRange)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    lossWindow = LossWindow(app)
    lossWindow.show()
    sys.exit(app.exec_())
