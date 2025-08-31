import os

from qtpy.QtCore import QSize
from qtpy.QtGui import QTextCursor, QDesktopServices
from qtpy.QtWidgets import QApplication, QSizePolicy, QWidget

from mfixgui.tools.qt import get_ui, get_icon

from mfixgui.version import get_version
from mfixgui.version_info import get_version_info


class AboutWidget(QWidget):
    """ Display version information for MFiX and its dependencies """

    def __init__(self):
        """ Create widget for About filemenu pane """
        super(AboutWidget, self).__init__()
        get_ui("about.ui", widget=self)
        self.label_about.setSizePolicy(
            QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum)
        )
        self.copy_btn.clicked.connect(self.copy_version_info)
        self.copy_btn.setIcon(get_icon('copy.svg'))

        self.textbrowser_version_info.setText(os.linesep.join(get_version_info()))

        size = self.label_supportforum.sizeHint().height()
        self.copy_btn.setIconSize(QSize(size, size))
        self.label_update_notification.setVisible(False)

    def copy_version_info(self):
        """ copy contents of TextBrowser to clipboard """
        txt = self.textbrowser_version_info.toPlainText()
        QApplication.instance().clipboard().setText(txt)
        self.select_text(txt)

    def select_text(self, txt):
        """ for visual feedback that the text was copied """
        cursor = self.textbrowser_version_info.textCursor()
        cursor.setPosition(0)
        cursor.setPosition(len(txt), QTextCursor.KeepAnchor)
        self.textbrowser_version_info.setTextCursor(cursor)

    def show_update_available(self, latest_version):
        self.label_update_notification.setText(
            f"You are running MFiX {get_version()}. "
            f"Version {latest_version} is now available for download from "
            f'<a href="https://mfix.netl.doe.gov/products/mfix/download/">the MFiX website</a>'
        )
        self.label_update_notification.setVisible(True)
