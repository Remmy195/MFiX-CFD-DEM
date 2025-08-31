"""Utility functions which depend on Qt"""

import copy
import os
import random
import sys

from qtpy import uic
from qtpy.QtCore import Qt, QSize, QSettings
from qtpy.QtGui import QColor, QIcon, QPixmap, QImage
from qtpy.QtWidgets import (QFrame,
                            QWidget,
                            QListWidget,
                            QStyleFactory,
                            QToolButton,
                            QStackedWidget)

from mfixgui.tools import SCRIPT_DIR

SETTINGS = QSettings("MFIX", "MFIX")

def set_item_noedit(item):
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)


def set_item_enabled(item, enabled):
    """Enable/disable items which do not have a setEnabled method, like menu
    items"""
    if not hasattr(item, 'flags'):
        return
    flags = item.flags()
    if enabled:
        flags |= Qt.ItemIsEnabled
    else:
        flags &= ~Qt.ItemIsEnabled
    item.setFlags(flags)


def item_enabled(item):
    flags = item.flags()
    return bool(flags & Qt.ItemIsEnabled)

def set_enabled(widget_or_item, enabled):
    # Do we need this?  Is it OK to just use set_item_enabled everywhere?
    if isinstance(widget_or_item, QWidget):
        widget_or_item.setEnabled(enabled)
    else:
        set_item_enabled(widget_or_item, enabled)

def get_widget_keys(w):
    keys = []
    if k := getattr(w, 'key', None):
        keys.append(k)
    if k:= getattr(w, 'keys', None):
        keys.extend(k)
    return keys

def get_combobox_item(cb, n):
    """Return the nth menu item from a combobox"""
    if not hasattr(cb, 'items'):
        cb.items = {}
    model = cb.model()
    item = model.item(n, 0)
    # Prevent Python objects from going out of scope, so we can attach attributes (e.g. 'tooltip0')
    cb.items[n] = item
    return item

def set_combobox_tooltip(combobox, index=None, default=''):
    """ Sets combobox tooltip to match item tooltip """
    if index is None:
        index = combobox.currentIndex()
    item = get_combobox_item(combobox, index)
    combobox.setToolTip(item.toolTip() if item else default)

def get_selected_row(table):
    """get index of selected row from a QTableWidget"""
    # note, currentRow can return  >0 even when there is no selection
    rows = set(i.row() for i in table.selectedIndexes())
    return None if not rows else rows.pop()


def get_selected_rows(table):
    """get index of selected row from a QTableWidget"""
    # note, currentRow can return  >0 even when there is no selection
    rows = set(i.row() for i in table.selectedIndexes())
    return sorted(list(rows))


def get_image_path(name):
    """"get path to images"""
    path = os.path.join(SCRIPT_DIR, 'images', name)
    # We do this because \ in stylesheet is problematic, see issues/1388
    if os.name == 'nt':
        path = path.replace('\\', '/')
    return path


PIXMAP_CACHE = {}
def get_pixmap(name, width, height):
    pixmap = PIXMAP_CACHE.get((name, width, height))
    if pixmap is None:
        pixmap = QPixmap(get_image_path(name)).scaled(
            width, height, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        PIXMAP_CACHE[(name, width, height)] = pixmap
    return pixmap


icon_cache = {}
def get_icon(name, resample=True):
    """Return image inside a QIcon object
    default: default image name or icon
    resample: if True, resample icon pixmaps for usual sizes
    (16, 24, 32, 48, 96, 128, 256). This is recommended for QMainWindow icons
    created from SVG images on non-Windows platforms due to a Qt bug (see
    http://code.google.com/p/spyderlib/issues/detail?id=1314)."""

    icon = icon_cache.get(name)
    if icon:
        return icon

    fname = get_image_path(name)
    if not os.path.exists(fname):
        for ext in '.svg', '.png':
            if os.path.exists(fname + ext):
                fname = fname + ext
                break

    icon = QIcon(fname)

    if resample:
        icon0 = QIcon()
        for size in (16, 24, 32, 48, 96, 128, 256, 512):
            icon0.addPixmap(icon.pixmap(size, size))
        icon = icon0
    icon_cache[name] = icon
    return icon


def widget_iter(widget):
    """iterator to recursively iterate over widgets"""
    if isinstance(widget, QStackedWidget):
        for n in range(widget.count()):
            yield from widget_iter(widget.widget(n))
        return
    for child in widget.children():
        if child.children():
            yield from widget_iter(child)
        yield child


def set_alignment(obj, align):
    # Some versions of Qt wrappers want int, some want QtAlignment
    try:
        obj.setAlignment(int(align))
    except TypeError:
        obj.setAlignment(align)

def random_pastel_color():
    """generate and return a random pastel color"""
    return [(random.randint(0, 128) + 100) / 255.0 for i in range(3)]


class CellColor:
    """
    A class to store color information and return '' if str or print is called
    on it. This is used to store colors in cells of a table.
    """

    def __init__(self, color=(1, 0, 0), text=''):

        if isinstance(color, (list, tuple)):
            self.color = QColor(*color)
        else:
            self.color = QColor(color)
        self.text = text

    @property
    def color_int(self):
        return [255 * c for c in self.color_float]

    @property
    def color_float(self):
        return self.color.getRgbF()[:3]

    @property
    def color_hex(self):
        return self.color.name()

    @property
    def qcolor(self):
        return self.color

    def __repr__(self):
        return self.text

    def __deepcopy__(self, memo):
        return CellColor(
            copy.deepcopy(self.color_hex), copy.deepcopy(self.text))

    def rand(self):
        self.color.setRgbF(*random_pastel_color())


def deepcopy_dict(dirty_dict, qobjects=False):
    '''deep copy a dictionary that has Qt objects in it
    Note: python 3.6+ can't copy Qt objects like QPixmap
    setting qobjects=True will copy the qt objects'''

    clean_dict = {}

    for key, value in dirty_dict.items():
        if isinstance(value, dict):
            clean_dict[key] = deepcopy_dict(value, qobjects)
        elif isinstance(value, QPixmap):
            if qobjects:
                clean_dict[key] = QPixmap(value)
        elif isinstance(value, QColor):
            if qobjects:
                clean_dict[key] = QColor(value)
        else:
            clean_dict[key] = copy.deepcopy(value)
    return clean_dict


def clear_layout(layout):
    """given a layout, clear all widgets"""
    while layout.count():
        item = layout.takeAt(0)
        if item is None:
            break
        widget = item.widget()
        if widget:
            widget.setParent(None)
            widget.deleteLater()


def get_separator(vertical=True, parent=None):
    """create a QFrame that looks like a separator"""
    f = QFrame
    line = f(parent)
    if vertical:
        line.setFrameShape(f.VLine)
    else:
        line.setFrameShape(f.HLine)
    line.setFrameShadow(f.Sunken)
    return line


_main_icon_size = QSize(24, 24)
_sub_icon_size = QSize(16, 16)


def set_main_icon_size(x):
    global _main_icon_size
    _main_icon_size = QSize(int(x), int(x))


def main_icon_size():
    return _main_icon_size


def set_sub_icon_size(x):
    global _sub_icon_size
    _sub_icon_size = QSize(int(x), int(x))


def sub_icon_size():
    return _sub_icon_size


def main_icon_height():
    return main_icon_size().height()


def sub_icon_height():
    return sub_icon_size().height()

def get_ui(ui_file, widget=None):
    """ returns the qt widget loaded from .ui file """
    ui_path = os.path.join(SCRIPT_DIR, 'uifiles', ui_file)
    if widget is None:
        return uic.loadUi(ui_path)
    return uic.loadUi(ui_path, widget)

def get_thumbnail(project_dir_path):
    thumb_path = os.path.join(project_dir_path, '.thumbnail')
    if os.path.isfile(thumb_path):
        pix = QPixmap()
        pix.convertFromImage(QImage(thumb_path, 'PNG'))
        icon = QIcon(pix)
    else:
        icon = get_icon('missing_thumbnail.png')
    return icon

def view_mode(is_tile_mode):
    return QListWidget.IconMode if is_tile_mode else QListWidget.ListMode


def get_preferred_style(cmdline_style, current_style):
    """ return valid style from cmdline, settings, or platform default """
    styles_avail = [s.lower() for s in QStyleFactory.keys()]

    if cmdline_style:
        if cmdline_style.lower() in styles_avail:
            return cmdline_style.lower()
        print("Style '%s' not available, ignoring" % cmdline_style)

    style_cf = SETTINGS.value("app_style")
    if style_cf:
        if style_cf.lower() in styles_avail:
            return style_cf.lower()
        print("Style '%s' not available, clearing from settings" % style_cf)
        try:
            SETTINGS.remove("app_style")
        except Exception as e:
            pass

    if "darwin" in sys.platform:
        return "macintosh"

    if "windows" in current_style and "fusion" in styles_avail:
        return "fusion"

    return current_style


def find_gui(parent):
    '''return the main gui instance by recursively searching up the inheritance
    tree'''
    # This is gross, gui should be available in global namespace
    if not hasattr(parent, 'objectName'):
        return None
    if parent.objectName() == 'mfixgui':
        return parent
    return find_gui(parent.parent())

def make_toolbutton(icon, callback, tooltip='', size=None):
    tb = QToolButton()
    tb.clicked.connect(callback)
    tb.setIcon(get_icon(icon))
    if size is not None:
        tb.setIconSize(size)
    tb.setToolTip(tooltip)
    return tb

def trim_table_height(tw):
    header_height = tw.horizontalHeader().height()
    scrollbar_height = tw.horizontalScrollBar().isVisible() * (4+tw.horizontalScrollBar().height())
    nrows = tw.rowCount()
    if nrows==0:
        row_height = 0
        height = header_height+scrollbar_height
    else:
        row_height = tw.rowHeight(0)
        height =  (header_height+scrollbar_height
                   + nrows*row_height + 4) # extra to avoid unneeded scrollbar
    tw.setMaximumHeight(height)
    # issues/1631
    tw.setMinimumHeight(header_height+ min(10, nrows)*row_height)
    tw.updateGeometry() #? needed?

def append_stylesheet(widget, ss):
    # XXX We're putting the new stylesheet entry at the start
    # rather than the end because otherwise font size of filemenu doesn't change
    if not ss.endswith(";"):
        ss += ';'
    term = ss.split(':')[0]
    # remove redundant settings
    s0 = widget.styleSheet().strip()
    if s0:
        s0 = '\n'.join(l for l in s0.splitlines() if not l.startswith(term))
    if s0:
        widget.setStyleSheet(ss + '\n' + s0)
    else:
        widget.setStyleSheet(ss)

def get_font_size():
    if int(SETTINGS.value('use_system_font', 0)):
        return None
    return int(SETTINGS.value('font_size', 12))
