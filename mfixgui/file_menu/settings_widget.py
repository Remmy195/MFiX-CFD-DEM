from PyQt5.QtWidgets import QSplitter, QStyleFactory, QWidget, QColorDialog, QToolButton
from PyQt5.QtCore import QSize
from PyQt5.QtGui import QColor, QFont

from mfixgui.animations import animation_speed
from mfixgui.tools.qt import get_ui, SETTINGS, widget_iter, append_stylesheet

class SettingsWidget(QWidget):
    def __init__(self, gui):
        """ Create widget for Settings filemenu pane """
        super(SettingsWidget, self).__init__()
        self._init_done = False
        get_ui("settings.ui", widget=self)
        self.gui = gui

        self.checkbox_developer_mode.clicked.connect(gui.enable_developer_mode)
        self.checkbox_collapse_splitters.clicked.connect(self.collapse_splitters)
        self.checkbox_show_resources.clicked.connect(gui.resource_monitor.show)
        self.groupbox_animation_speed.clicked.connect(lambda n: self.spinbox_animation_speed.setValue(200 if n else 0))
        self.spinbox_animation_speed.valueChanged.connect(lambda n: SETTINGS.setValue('animation_speed', n))
        self.checkbox_screenshot_res.clicked.connect(lambda state: SETTINGS.setValue('enable_screenshot_res', int(state)))
        self.combobox_output_template.currentIndexChanged.connect(self.set_template)
        self.combobox_sms.activated.connect(gui.set_sms_enabled)
        self.label_vtk_color1.mousePressEvent = lambda ev: self.change_color(1)
        self.label_vtk_color2.mousePressEvent = lambda ev: self.change_color(2)
        self.spinbox_font_size.setKeyboardTracking(False)
        self.spinbox_font_size.valueChanged.connect(self.change_font_size)
        self.checkbox_use_system_font.clicked.connect(self.set_use_system_font)

        # GUI style
        current_style = gui.style
        avail_styles = [s.lower() for s in QStyleFactory.keys()]
        if current_style and current_style not in avail_styles:
            gui.error("Application style %s not available" % current_style, popup=True)
            avail_styles.append(current_style) # ?

        cb = self.combobox_appstyle
        cb.clear()
        cb.addItems(avail_styles)
        cb.currentIndexChanged.connect(lambda n, cb=cb: SETTINGS.setValue('app_style', cb.currentText()))

        self.pushbutton_reset.clicked.connect(self.reset_settings)
        self.update()
        self._init_done = True


    def reset_settings(self):
        ok = self.gui.message(text='This will reset all MFiX global settings\n'
                                   'to their original default values.\n'
                                   'Are you sure?',
                              buttons=['ok','cancel'],
                              default='cancel')
        if ok == 'ok':
            self.gui.print_internal("Settings cleared.", color='blue')
            SETTINGS.clear()
            SETTINGS.sync()
            self.update()
            # How many other things do we need to reset here?
            self.gui.open_widget.clear_recent()
            self.gui.set_sms_enabled(False) # Default

    def change_font_size(self, n):
        self.font_size = n
        use_sys = int(SETTINGS.value("use_system_font", 0))
        if use_sys or n is None:
            # Clear out stylesheet
            for x in self.gui.ui, self.gui.file_menu:
                ss = x.styleSheet()
                x.setStyleSheet('\n'.join(l for l in ss.splitlines()
                                          if not 'font-size' in l))
            self.gui.change_font('', None)
        else:
            self.gui.change_font('', n)
            append_stylesheet(self.gui.ui,("font-size: %spt" % n))
            append_stylesheet(self.gui.file_menu, ("font-size: %spt;" % n))
        self.gui.set_file_menu_width()
        if n:
            SETTINGS.setValue('font_size',n)
            SETTINGS.sync()

    def set_use_system_font(self, b):
        SETTINGS.setValue('use_system_font', int(b))
        SETTINGS.sync()

        if b:
            # Clear out stylesheet
            for x in self.gui.ui, self.gui.file_menu:
                ss = x.styleSheet()
                x.setStyleSheet('\n'.join(l for l in ss.splitlines()
                                          if not 'font-size' in l))
            self.change_font_size(None)
        else:
            font_size = SETTINGS.value('font_size', 12)
            if not font_size:
                font_size = 12
            font_size = int(font_size)
            self.change_font_size(font_size)

        self.spinbox_font_size.setEnabled(not b)
        self.label_font_size.setEnabled(not b)


    def set_template(self):
        template_name = self.combobox_output_template.currentText()
        SETTINGS.setValue('template', template_name)
        self.gui.set_unsaved_flag() # Allow re-save with new template
        if template_name == 'Custom' and self._init_done:
            self.gui.warning("Custom template must be supplied by user.  See MFiX documentation.",
                             popup=True)

    def update(self):
        self.checkbox_developer_mode.setChecked(
            int(SETTINGS.value('developer_mode', 0)))
        self.checkbox_collapse_splitters.setChecked(
            int(SETTINGS.value('collapse_qsplitter', 0)))
        self.checkbox_show_resources.setChecked(
            int(SETTINGS.value('show_resources', 0)))
        speed = animation_speed()
        self.groupbox_animation_speed.setChecked(speed > 0)
        self.spinbox_animation_speed.setValue(speed)
        self.checkbox_screenshot_res.setChecked(
            int(SETTINGS.value('enable_screenshot_res', 0)))
        font_size = SETTINGS.value('font_size', 12) or 12
        font_size = int(font_size)
        self.spinbox_font_size.setValue(font_size)

        use_sys = int(SETTINGS.value('use_system_font', 0))
        self.checkbox_use_system_font.setChecked(use_sys)
        self.label_font_size.setEnabled(not use_sys)
        self.spinbox_font_size.setEnabled(not use_sys)
        self.label_font_size.hide()

        # MFX file template
        current_template = SETTINGS.value('template', 'Standard')
        cb = self.combobox_output_template
        cb.setCurrentText(current_template)

        # SMS mode
        sms = SETTINGS.value('SMS', 0)
        if sms in ('False', 'false'):
            sms = 0
        elif sms in ('True', 'true'):
            sms = 1
        sms = int(sms)
        self.combobox_sms.setCurrentIndex(1 if sms else 0)
        self.combobox_appstyle.setCurrentText(self.gui.style)

        # VTK colors
        self.set_color(self.label_vtk_color1,
                              QColor(SETTINGS.value('vtk_bck_color1', '#aaffff')))

        self.set_color(self.label_vtk_color2,
                              QColor(SETTINGS.value('vtk_bck_color2', '#ffffff')))


    def collapse_splitters(self, enable):
        for widget in widget_iter(self.gui):
            if isinstance(widget, QSplitter):
                widget.setChildrenCollapsible(enable)

        SETTINGS.setValue("collapse_qsplitter", int(enable))

    def change_color(self, ncolor):
        key = "vtk_bck_color%s"%ncolor
        initial = QColor(SETTINGS.value(key, '#aaffff' if ncolor==1 else '#ffffff'))
        color = QColorDialog.getColor(
            initial=initial, parent=self, title='Select background color')

        label = self.label_vtk_color1 if ncolor==1 else self.label_vtk_color2
        if color.isValid():
            self.set_color(label, color)
            SETTINGS.setValue(key, color.name())

            c1 = QColor(SETTINGS.value('vtk_bck_color1', '#aaffff'))
            c2 = QColor(SETTINGS.value('vtk_bck_color2', '#ffffff'))
            #push update
            for tab in self.gui.graphics_tabs:
                if tab.vtk:
                    tab.vtk_widget.set_background_color(c1.getRgbF()[:3], c2.getRgbF()[:3])
            self.gui.vtkwidget.set_background_color(c1.getRgbF()[:3], c2.getRgbF()[:3])

    def set_color(self, widget, color):
        rgb = tuple(int(255 * c) for c in color.getRgbF()[:3])
        ss = "background-color: rgb(%s,%s,%s); border: 1px solid black; border-radius: 2px;"%rgb
        widget.setStyleSheet(ss)
