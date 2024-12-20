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

class Settings:

    def __init__(self):
        pass

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

    