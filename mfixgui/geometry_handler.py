# Methods for geometry task pane window

from mfixgui.constants import *
from mfixgui.tools.qt import widget_iter
from mfixgui.widgets.base import BaseWidget

class GeometryHandler:
    def init_geometry(self):
        ui = self.ui.geometry
        for w in widget_iter(ui):
            if isinstance(w, BaseWidget):
                w.post_update = self.clear_mesh_accepted
        ui.checkbox_keyword_no_k.post_update = self.handle_no_k
        ui.pushbutton_mesh_autosize.keys = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']
        ui.pushbutton_mesh_autosize.lock_button = True # see gui.enable_input

    def enable_z_input(self, enable):
        ui = self.ui.geometry
        for w in (ui.lineedit_keyword_z_min,
                  ui.lineedit_keyword_z_max,
                  ui.label_Z):
            self.set_widget_enabled(w, enable)

    def reset_geometry(self):
        self.enable_z_input(True)

    def setup_geometry(self, allow_disabled_tab=False):
        ui = self.ui.geometry
        no_k = self.project.get_value('no_k', default=False)
        # Issues/871
        cb = ui.checkbox_keyword_no_k
        if self.project.solver == PIC:
            reason = "PIC model does not support 2D simulation"
            if no_k:
                self.error(reason, popup=True)
                self.unset_keyword('no_k')
            self.set_widget_enabled(cb, False, reason)
            cb.setChecked(False)
            self.enable_z_input(True)
        else:
            self.set_widget_enabled(cb, True)
            cb.setChecked(no_k)
            disable = no_k and self.project.solver in (DEM,CGP,SQP,GSP)
            self.enable_z_input(not disable)


    def set_z_max_from_d_p0(self):
        mmax = self.project.get_value('mmax', default=1)
        d =[self.project.get_value('d_p0', default=0, args=[i])
            for i in range(1, 1+mmax)]
        if d:
            self.update_keyword('z_min', 0)
            self.update_keyword('z_max', max(d))


    def handle_no_k(self):
        self.clear_mesh_accepted()
        no_k = self.project.get_value('no_k', default=False)
        if no_k:
            self.retain_keyword('kmax')
            self.unset_keyword('kmax')
        else:
            val = self.get_retained_keyword('kmax')
            if val is not None:
                self.update_keyword('kmax', val)
        if no_k and self.project.solver in (DEM,CGP,SQP,GSP):
            self.enable_z_input(False)
            self.set_z_max_from_d_p0()
        else:
            self.enable_z_input(True)
        self.update_background_mesh()
