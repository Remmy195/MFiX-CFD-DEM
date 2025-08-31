import os

from qtpy.QtCore import Qt
from qtpy.QtGui import QDesktopServices, QIcon

from qtpy.QtWidgets import QLabel, QListWidgetItem, QWidget

from mfixgui.tools import SCRIPT_DIR, get_mfix_doc_html
from mfixgui.tools.qt import get_ui

from mfixgui.version_info import MFIX_VERSION

try:
    from nodeworks.tools import get_nodeworks_doc
    from nodeworks import __version__ as NODEWORKS_VERSION
    NODEWORKS_DOC_PATH = get_nodeworks_doc()
except ImportError:
    NODEWORKS_VERSION = NODEWORKS_DOC_PATH = None


class HelpWidget(QWidget):
    def __init__(self, mfixgui):
        """ Create widget for help filemenu pane """
        super(HelpWidget, self).__init__()
        get_ui("help.ui", widget=self)
        self.gui = mfixgui

        docs = get_mfix_doc_html()
        mfix_docs = ("#" if docs is None
                     else "file://" + docs + "/index.html")

        self.label_mfix_local.setText(
            f'<a href="{mfix_docs}">' f"Local MFiX Documentation (v {MFIX_VERSION})</a>")

        if NODEWORKS_VERSION:
            if NODEWORKS_DOC_PATH:
                self.label_nodeworks_local.setText('<a href="file://%s"">Local Nodeworks Documentation (v. %s})</a>' %
                                                   (NODEWORKS_DOC_PATH, NODEWORKS_VERSION))
                self.label_nodeworks_local.setVisible(True)
            else:
                self.label.nodeworks_local.setVisible(False)
            self.label_nodeworks_online.setVisible(True)
            self.label_nodeworks_forum.setVisible(True)
        else:
            self.label_nodeworks_local.setVisible(False)
            self.label_nodeworks_online.setVisible(False)
            self.label_nodeworks_forum.setVisible(False)

        lw = self.listwidget_tutorial
        lw.setAttribute(Qt.WA_MacShowFocusRect, 0)
        lw.setFrameStyle(lw.NoFrame)

        for thumb, text in TUTORIALS:
            label = QLabel(text)
            label.setOpenExternalLinks(True)
            item = QListWidgetItem()
            item.setIcon(QIcon(thumb))
            item.setFlags(Qt.ItemIsEnabled)
            lw.addItem(item)
            lw.setItemWidget(item, label)


def make_tutorial(tut):
    image, title, tutorial, youtube, _ = tut
    thumb = os.path.join(SCRIPT_DIR, "images", image)
    docs = get_mfix_doc_html()
    url = ("file://"+docs+"/tutorials/"+tutorial if docs
           else "#")
    video = "Video" if youtube is None else f'<a href="{youtube}">Video</a>'
    text = f'<b>{title}</b><p><a href="{url}">Text</a> | {video}</p>'
    return (thumb, text)


# URLs and thumbnails for Tutorials in sphinx documentation
TUTORIALS = map(
    make_tutorial,
    [
        (
            "gui_tfm_2d_thumbnail.png",
            "Two-dimensional fluid bed two-fluid model (TFM)",
            "tutorial_tfm.html",
            "https://youtu.be/rZgdGH2pkx4",
            ["2d", "tfm"],
        ),
        (
            "gui_dem_2d_thumbnail.png",
            "Two-dimensional fluid bed discrete-element model (DEM)",
            "tutorial_dem.html",
            None,
            ["2d", "dem"],
        ),
        (
            "gui_sphere_thumbnail.png",
            "3D Single-phase flow over a sphere",
            "tutorial_sphere.html",
            None,
            ["3d", "single", "cartesian"],
        ),
        (
            "gui_3d_fluidbed_thumbnail.png",
            "3D Fluid Bed (TFM & DEM)",
            "tutorial_3d_fluidbed.html",
            None,
            ["3d", "tfm", "cartesian"],
        ),
        (
            "gui_hopper_thumbnail.png",
            "3D Hopper discrete-element model (DEM)",
            "tutorial_hopper.html",
            "https://youtu.be/kEiPs-phndg",
            ["3d", "dem", "cartesian"],
        ),
        (
            "gui_dem_chutes_thumbnail.png",
            "DEM Granular Flow Chutes",
            "tutorial_dem_chutes.html",
            None,
            ["3d", "dem", "cartesian"],
        ),
    ],
)
