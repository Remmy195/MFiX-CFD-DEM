#!/usr/bin/env python
# Dialog for ramp_* keys.
# Only implemented for BC_* keys, for now

import sys

from qtpy.QtWidgets import QApplication, QDialog

from mfixgui.widgets.base import LineEdit
from mfixgui.tools.qt import get_ui
from mfixgui.tools.keyword_args import mkargs

class RampPopup(QDialog):
    def __init__(self, parent=None):
        super(RampPopup, self).__init__(parent)
        self.parent = parent
        ui = self.ui = get_ui('ramp_popup.ui', self)
        if parent:
            LineEdit.report_value_error = parent.popup_value_error
            LineEdit.report_value_required = parent.popup_value_required

        for le in (ui.lineedit_t0, ui.lineedit_t1):
            le.min = 0.0

        for le in (ui.lineedit_t0, ui.lineedit_t1,
                   ui.lineedit_v0, ui.lineedit_v1):
            le.dtype = float
            le.allow_parameters = True
            le.value_updated.connect(self.update)
        ui.button_reset.clicked.connect(self.reset)
        ui.button_ok.clicked.connect(self.close)


    def reset(self):
        ui = self.ui
        ui.lineedit_t0.setText('0.0')
        ui.lineedit_t1.setText('')
        ui.lineedit_v0.setText('')
        ui.lineedit_v1.setText('')
        for suff in ('t0','t1','v0','v1'):
            ramp_key = self.key + '_ramp_' + suff
            if self.parent:
                gui = self.parent
                for BC in self.indices:
                    args = mkargs(ramp_key, bc=BC, phase=self.phase)
                    gui.retain_keyword(ramp_key, args=args)
                    gui.unset_keyword(ramp_key, args=args)


    def update(self, *args):
        if not self.parent:
            print(args)
            return
        gui = self.parent
        w, value_dict, _ = args
        key, val = value_dict.popitem()
        for BC in self.indices:
            args = mkargs(key, bc=BC, phase=self.phase)
            gui.update_keyword(key, val, args=args)
            if key.endswith('_ramp_v0'):
                gui.update_keyword(key[:-8], val, args=args)

    def popup(self, key=None, unit='', indices=[], phase=None):
        #print("POPUP", key, indices)

        if not indices: # nothing to act on
            return
        self.key = key
        self.indices = indices
        self.phase = phase
        gui = self.parent
        proj = gui.project
        ui = self.ui
        ui.label_unit0.setText(unit)
        ui.label_unit1.setText(unit)

        BC0 = indices[0]
        for suff in ('t0','t1','v0','v1'):
            ramp_key = key + '_ramp_' + suff
            le = getattr(ui, 'lineedit_' + suff)
            le.key = ramp_key
            gui.add_tooltip(le, key=ramp_key)
            args = mkargs(ramp_key, bc=BC0, phase=phase)
            val = proj.get_value(ramp_key, args=args)
            #print("GET", ramp_key, args, 'RETURNS', val)
            default = proj.keyword_doc.get(key, {}).get('initpython')
            le.default(default)
            #print(key, "DEFAULT", default)
            if val is None:
                val = None
                if ramp_key.endswith('_v0'):
                    val = proj.get_value(key, args=args)
                    if val is None:
                        val = gui.get_retained_keyword(key, args=args,
                                                       default=0.0 if suff=='t0' else None)
                if val is None:
                    val = gui.get_retained_keyword(ramp_key, args=args,
                                                   default=0.0 if suff=='t0' else None)
                if val is not None:
                    for BC in indices:
                        gui.update_keyword(ramp_key, val,
                                           args=mkargs(ramp_key, bc=BC, phase=phase))
            if val is None:
                le.setText('' if default is None else str(default))
            else:
                le.setText(str(val))


        ui.lineedit_t0.setFocus()
        self.show()
        self.raise_()
        self.activateWindow()

def main():
    args = sys.argv
    qapp = QApplication(args)
    ramp_popup = RampPopup()
    ramp_popup.show()

    qapp.exec_()
    qapp.deleteLater()

    sys.exit()

if __name__ == '__main__':
    main()
