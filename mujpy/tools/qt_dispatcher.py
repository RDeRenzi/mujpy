from PyQt5 import QtCore, QtWidgets
import traceback
import io
from contextlib import redirect_stdout, redirect_stderr

class QtDispatcher(QtCore.QObject):
    """Run mudashed plotting operations on Qt's GUI thread."""

    requested = QtCore.pyqtSignal(object)

    def __init__(self, log=None):
        super().__init__()
        self._log = log

        app = QtWidgets.QApplication.instance()
        if app is None or self.thread() != app.thread():
            raise RuntimeError(
                "Create QtDispatcher on the GUI thread, after %matplotlib qt"
            )

        self.requested.connect(
            self._execute,
            QtCore.Qt.QueuedConnection,
        )

#    @QtCore.pyqtSlot(object)
#    def _execute(self, callback):
#        try:
#            callback()
#        except Exception:
#            message = traceback.format_exc()
#            if self._log is not None:
#                self._log(message)
#            else:
#                traceback.print_exc()

    def submit(self, callback):
        self.requested.emit(callback)




    @QtCore.pyqtSlot(object)
    def _execute(self, callback):
        if self._log is None:
    
            try:
                callback()
            except Exception:
                traceback.print_exc()
            return
    
        output = io.StringIO()
        try:
            with redirect_stdout(output), redirect_stderr(output):
                try:
                    callback()
                except Exception:
                    traceback.print_exc()
        finally:
            text = output.getvalue()
            if text:
                self._log(text)

