import pyqtgraph as pg

import sys

import numpy as np
from scipy.fft import fft

from superqt import QLabeledRangeSlider, QLabeledSlider, QDoubleSlider

import matplotlib.transforms as transforms
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt, QEvent

import loss_window
import window
import transmit_window

from scipy import signal

if __name__ == "__main__":
    qapp = QtWidgets.QApplication.instance()

    if not qapp:
        qapp = QtWidgets.QApplication(sys.argv)

    app = window.ApplicationWindow()
    plot = transmit_window.TransmitWindow()
    loss = loss_window.LossWindow()

    app.show()
    plot.show()

    # loss.show()

    app.activateWindow()
    app.raise_()
    qapp.exec()
