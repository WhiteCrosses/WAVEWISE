from superqt import QLabeledRangeSlider, QLabeledSlider  # type: ignore
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvas
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt  # type: ignore


class LossWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.createWidgets()

    def createWidgets(self):
        self.mainLayout = QtWidgets.QFormLayout()

        self.rangeSelectorSlider = QLabeledRangeSlider(Qt.Horizontal)
        self.rangeSelectorSlider.setMinimum(70)
        self.rangeSelectorSlider.setMaximum(6000)
        self.rangeSelectorSlider.setSingleStep(1)
        self.rangeSelectorSlider.valueChanged.connect(self.rangeChangeSlider)

        self.startFreqBox = QtWidgets.QDoubleSpinBox()
        self.startFreqBox.setMinimum(70)
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
        # self.runButton.clicked.connect(self.run)

        self.startFreqUnitBox = QtWidgets.QComboBox()
        self.endFreqUnitBox = QtWidgets.QComboBox()

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)

        self.mainLayout.addRow(self.canvas)
        self.mainLayout.addRow("Range to scan [MHz]", self.rangeSelectorBox)
        self.mainLayout.addRow("Single step size [MHz]", self.stepSelector)
        self.mainLayout.addRow("Select gain [mdB]", self.gainSelector)
        self.mainLayout.addRow(self.runButton)

        self.setLayout(self.mainLayout)

    def rangeChangeSlider(self):
        self.startFreqBox.setValue(self.rangeSelectorSlider.value()[0])
        self.endFreqBox.setValue(self.rangeSelectorSlider.value()[1])
        self.selectedRange = self.rangeSelectorSlider.value()

    def rangeChangeBox(self):
        self.selectedRange = (self.startFreqBox.value(),
                              self.endFreqBox.value())
        self.rangeSelectorSlider.setValue(self.selectedRange)
