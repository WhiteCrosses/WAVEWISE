from superqt import QLabeledRangeSlider, QLabeledSlider  # type: ignore
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvas
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from PyQt5.QtCore import Qt  # type: ignore
import numpy as np
import threading
import time
from scipy import signal
from scipy.fft import fftshift, fftfreq
from scipy import fft
import adi  # type: ignore

# TODO detecting sent signal works. it prints value of received signal. now iteration through range of frequencies and displaying them

class NoSDR:
    def __init__(self):
        self.sig = False
        self.rx_data = np.zeros_like(np.shape(5000))

    def rx(self):
        if self.sig:
            return self.rx_data
        else:
            return np.zeros_like(self.rx_data)


    def tx(self, signal):
        self.sig = True
        self.rx_data = signal

    def tx_destroy_buffer(self):
        self.sig = False


class LossWindow(QtWidgets.QWidget):
    def __init__(self, app):
        super().__init__()
        self.app = app
        self.selectedRange = (200,3000)
        self.sampleRate = int(10e6)
        self.n_samples = 5000
        self.center_freq = int(100e6)
        self.bufferSize = 5000

        self.plutoInit()
        self.result = np.array([])
        self.currFreq = 0
        self.createWidgets()

    def plutoInit(self):
        gain = 50.0  # allowable range is 0 to 74.5 dB
        self.sdr = NoSDR()

        self.sdr.sample_rate = self.sampleRate
        
        self.sdr.rx_rf_bandwidth = self.sampleRate
        self.sdr.rx_lo = self.center_freq
        self.sdr.rx_buffer_size = self.bufferSize
        self.sdr.rx_rf_bandwidth = int(self.sampleRate)
        self.sdr.rx_hardwaregain_chan0 = gain  # set receive gain
        
        self.sdr.tx_rf_badnwidth = self.sampleRate
        self.sdr.tx_hardwaregain_chan0 = -20
        self.sdr.gain_control_mode_chan0 = "manual"  # turn off AGC

    def createWidgets(self):
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
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def normalize(self, arr):
        # Find the maximum absolute magnitude in the array
        max_magnitude = np.max(np.abs(arr))
        # Divide the array by the maximum magnitude
        normalized_arr = arr / max_magnitude
        return normalized_arr

    def transmit(self, selectedFreq, gain):
        
        Fs = self.sampleRate

        t = np.arange(self.n_samples)/self.sampleRate
        
        selectedFreq = int(selectedFreq * 1e6)
        centerFreq = self.center_freq

        print("\n=====================================")
        print("selected freq: " + str(selectedFreq) + " Hz")

        self.samples = 0.5*np.exp(2.0j*np.pi*1e6*t)
        self.samples += 0.5*np.exp(2.0j*np.pi*2e6*t)
        self.samples += 0.5*np.exp(2.0j*np.pi*3e6*t)
        self.samples += 0.5*np.exp(2.0j*np.pi*4e6*t)
        self.samples = self.normalize(self.samples)
        self.samples *= 2**14

        self.sdr.rx_lo = int(centerFreq)        
        self.sdr.tx_lo = int(centerFreq)

        self.freqs_array = np.fft.fftshift(np.fft.fftfreq(self.n_samples, d=1/self.sampleRate)+centerFreq)
        
        sigOff = np.array([])
        sigOn = np.array([])

        for n in range (0, 1):
            for x in range(0, 5):
                raw_data = self.sdr.rx()
            off_raw = self.sdr.rx()

            off_fft = (np.fft.fft(off_raw))
            off_fft_db = 20*np.log10((np.abs(off_fft/self.n_samples)*2.0))
            
            
            off_fft_db = np.fft.fftshift(off_fft_db)
            nearest_f = self.find_nearest(self.freqs_array, selectedFreq)
            nearest_idx = np.where(self.freqs_array == nearest_f)[0][0]
            sigOff = (off_fft_db[nearest_idx])


        self.sdr.tx_cyclic_buffer = True
        self.sdr.tx(self.samples)
        for n in range (0, 1):
            for x in range(0, 5):
                raw_data = self.sdr.rx()
            on_raw = self.sdr.rx()
            on_fft = (np.fft.fft(on_raw))
            on_fft_db = 20*np.log10((np.abs(on_fft/self.n_samples)*2.0))
            on_fft_db = np.fft.fftshift(on_fft_db)

            on_peaks = signal.find_peaks(on_fft_db, height=-20,distance=100)
            nearest_peak = self.find_nearest(on_peaks[0], selectedFreq)
            nearest_peak_freq = self.find_nearest(self.freqs_array, nearest_peak)
            nearest_peak_freq_idx = int(0.7*self.n_samples) #np.where(self.freqs_array == nearest_peak_freq)[0][0]

            sigOn = (on_fft_db[nearest_peak_freq_idx])
        

        np.abs(on_fft/self.n_samples)
        diff = np.abs(np.fft.fftshift(on_fft)/self.n_samples) - np.abs(np.fft.fftshift(off_fft)/self.n_samples)
        loss = 20*np.log10((diff[nearest_peak_freq_idx]*2.0))

        self.sdr.tx_destroy_buffer()
        self.sdr.tx_cyclic_buffer = False

        print("SigOn: " + str(sigOn))
        print("SigOff: " + str(sigOff))
        print("idx: " + str(nearest_idx))
        print("idx2: "+ str(nearest_peak_freq_idx))

        #ret_on = 20*np.log10((np.abs(np.average(sigOn)/self.n_samples)*2.0))
        #ret_off = 20*np.log10((np.abs(np.average(sigOff)/self.n_samples)*2.0))

        return sigOn, sigOff, loss, on_fft_db, off_fft_db
        # append to self.result value of freq selected

    # TODO
    # wait for received fft, check peak values

    def run(self):
        self.currFreq = self.selectedRange[0]

        window = np.hanning(5)

        print(window)

        freq_array = np.arange(
            self.selectedRange[0]*1e6, self.selectedRange[1]*1e6, int(self.stepSelector.value()*1e6))

        idx = 0
        self.result = np.array([])
        self.result2 = np.array([])
        self.result3 = np.array([])

        # set click listener
        self.figure.canvas.mpl_connect('button_press_event', self.onclick)

        for i in range(self.selectedRange[0], self.selectedRange[1], self.stepSelector.value()):
            self.center_freq = int(i*1e6 - 1e6)
            sigOn, sigOff, loss, on_fft_db, off_fft_db = self.transmit(i, self.gainSelector.value())

            self.result = np.append(self.result, loss)
            self.result2 = np.append(self.result2, sigOn)
            self.result3 = np.append(self.result3, sigOff)

            #displayed = np.convolve(self.result, window, mode='same')
            #displayed2 = np.convolve(self.result2, window, mode='same')
            #displayed3 = np.convolve(self.result3, window, mode='same')

            displayed = self.result#np.convolve(self.result, window, mode='same')
            displayed2 = self.result2#np.convolve(self.result2, window, mode='same')
            displayed3 = self.result3#np.convolve(self.result3, window, mode='same')

            idx += 1
            

            if (idx > 10):
                self.figure.clear()
                self.ax = self.figure.add_subplot(311)
                self.ax2 = self.figure.add_subplot(312)
                self.ax3 = self.figure.add_subplot(313)

                self.ax.plot(freq_array[0:idx], displayed, 'r')
                self.ax.plot(freq_array[0:idx], displayed2, 'g')
                self.ax.plot(freq_array[0:idx], displayed3, 'b')


                self.ax2.plot(off_fft_db)
                self.ax3.plot(on_fft_db)

                def fmt(x, pos): return '{:.0f}'.format((x)/1e6, pos)
                self.ax.xaxis.set_major_formatter(
                    ticker.FuncFormatter(fmt))

                self.ax.set_xlim(
                    ((self.selectedRange[0] - 5)*1e6, (self.selectedRange[1] + 5)*1e6))

                self.canvas.draw()
                self.canvas.flush_events()

    def onclick(self, event):
        if event.button == 'q':
            quit()

    def rangeChangeSlider(self):
        self.startFreqBox.setValue(self.rangeSelectorSlider.value()[0])
        self.endFreqBox.setValue(self.rangeSelectorSlider.value()[1])
        self.selectedRange = self.rangeSelectorSlider.value()

    def rangeChangeBox(self):
        self.selectedRange = (self.startFreqBox.value(),
                              self.endFreqBox.value())
        self.rangeSelectorSlider.setValue(self.selectedRange)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    lossWindow = LossWindow(app)
    lossWindow.show()
    sys.exit(app.exec_())
