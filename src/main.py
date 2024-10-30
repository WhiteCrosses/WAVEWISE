import sys
from matplotlib.backends.qt_compat import QtWidgets
import matplotlib.pyplot as plt

import loss_window
import window
import transmit_window


if __name__ == "__main__":
    qapp = QtWidgets.QApplication.instance()

    if not qapp:
        qapp = QtWidgets.QApplication(sys.argv)

    app = window.ApplicationWindow()
    plot = transmit_window.TransmitWindow(app)
    app.plot = plot
    loss = loss_window.LossWindow()

    app.show()
    plot.show()

    # loss.show()

    app.activateWindow()
    app.raise_()
    qapp.exec()
