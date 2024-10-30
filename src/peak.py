import pyqtgraph as pg

import sys
import adi

import numpy as np
from scipy.fft import fft
# Krzywe transmisyjne

from superqt import QLabeledRangeSlider, QLabeledSlider, QDoubleSlider

import matplotlib.transforms as transforms
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt, QEvent

import time
import threading
import os
import peak
import csv


from scipy import signal


class Peak:
    def __init__(self, frequency, power, distance):
        self.frequency = frequency
        self.distance = distance
        self.min = self.frequency - self.distance
        self.max = self.frequency + self.distance
        self.distanceToNext = None
        self.tickCounter = 0
        self.isChecked = False
        self.power = power
        self.color = 'red'

    def found(self, x, power):
        self.frequency = x
        self.power = power
        self.min = self.frequency - self.distance
        self.max = self.frequency + self.distance
        self.resetTickCounter()

    def resetTickCounter(self):
        self.isChecked = True
        self.tickCounter = 0
