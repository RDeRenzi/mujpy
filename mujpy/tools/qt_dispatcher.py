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

    def submit(self, callback):
        self.requested.emit(callback)

    @QtCore.pyqtSlot(object)
    def _execute(self, callback):
        if self._log is None:
    
            try:
                self._run_and_focus(callback)
            except Exception:
                traceback.print_exc()
            return
    
        output = io.StringIO()
        try:
            with redirect_stdout(output), redirect_stderr(output):
                try:
                    self._run_and_focus(callback)
                except Exception:
                    traceback.print_exc()
        finally:
            text = output.getvalue()
            if text:
                self._log(text)
                
    def _run_and_focus(self, callback):
        import shutil
        import subprocess
        import traceback

        fig = callback()
        canvas = getattr(fig, "canvas", None)
        manager = getattr(canvas, "manager", None)
        window = getattr(manager, "window", None)

        if window is None:
            return fig

        def focus():
            try:
                if window.isMinimized():
                    window.showNormal()
                else:
                    window.show()

                window.raise_()
                window.activateWindow()

                wmctrl = shutil.which("wmctrl")
                if wmctrl and QtWidgets.QApplication.platformName() == "xcb":
                    subprocess.run(
                        [wmctrl, "-i", "-a", hex(int(window.winId()))],
                        check=True,
                        timeout=2,
                    )

            except Exception:
                if self._log is not None:
                    self._log(traceback.format_exc())
                else:
                    traceback.print_exc()

        QtCore.QTimer.singleShot(100, focus)
        return fig
