# Methods to deal with solids pic tab, split off from solids_handler.py

# Particle in Cell Model Task Pane Window

from .constants import *

PIC_defaults = ( # key, val, min, max, exclude_min, exclude_max
    ('des_epg_clip', 0.42, 0.0, 1.0, True, True),                 # (0.0, 1.0)
    ('fric_exp_pic', 3.0, 1.0, 5.0, False, False),           # [1.0, 5.0]
    ('psfac_fric_pic', 100.0, 1.0, 10000.0, False, False),   # [1.0, 10000.0]
    ('mppic_coeff_en1', 0.85, 0.0, 1.0, True, False),        # (0.0, 1.0]
    ('fric_non_sing_fac', 1.0e-7, 0.0, 1.0e-4, True, False), # (0.0, 1.0e-4]
    ('mppic_coeff_en_wall', 0.85, 0.0, 1.0, True, False),    # (0.0, 1.0]
    ('mppic_coeff_et_wall', 1.0, 0.0, 1.0, True, False),     # (0.0, 1.0]
    ('mppic_velfac_coeff', 1.0, 0.0, 1.0, False, False))     # [0.0, 1.0]

class SolidsPIC:

    def init_solids_pic(self):
        ui = self.ui.solids

        # Set up ranges for PIC variables.
        for (key, default, min, max, exclude_min, exclude_max) in PIC_defaults:
            le = getattr(ui, 'lineedit_keyword_%s'%key)
            le.min = min
            le.max = max
            le.exclude_min = exclude_min
            le.exclude_max = exclude_max

        # EP_STAR is REQUIRED
        # 1538 ep_star replaced by des_epg_clip
        for le in (ui.lineedit_keyword_des_epg_clip,
                   ui.lineedit_keyword_ep_star):
            default = 0.42
            le.required = True
            le.saved_value = default
            le.initpython = default

        le = ui.lineedit_keyword_pic_cfl
        le.post_update = self.setup_solids_pic_tab
        cb = ui.combobox_pic_cfl_control
        key = 'pic_cfl_control'
        cb.activated.connect(lambda idx:
                             self.update_keyword('pic_cfl_control',
                                                 CFL_CONTROL_TYPES[idx]))
        self.add_tooltip(cb, key=key)
        cb = ui.checkbox_keyword_pic_collision_damping
        cb.post_update = self.handle_pic_collision_damping

    def setup_solids_pic_tab(self):
        ui = self.ui.solids
        # We might be just visiting due to locate_keyword
        ui.PIC.setEnabled(ui.pushbutton_solids_PIC.isEnabled())
        # PSD table is shared between DEM and PIC per #1356 (should it be in Materials?)
        gb = ui.groupbox_psd
        par = gb.parent()
        if par != ui.frame_pic: # Widget has moved to DEM pane
            gb.hide()
            par.layout().removeWidget(gb)
            layout = ui.frame_pic.layout()
            layout.insertWidget(2, gb)
            gb.show()
        self.update_solids_psd_table()
        self.fixup_solids_table(ui.tablewidget_psd)

        key = 'pic_cfl_control'
        val = self.project.get_value(key, default='MAX')
        if val not in CFL_CONTROL_TYPES:
            self.error("Invalid value for key %s:  '%s', setting to MAX" %
                       (key, val),
                       popup=True)
            val = 'MAX'
            self.update_keyword(key, val)
        cb = ui.combobox_pic_cfl_control
        cb.setCurrentIndex(CFL_CONTROL_TYPES.index(val))
        enable =  self.project.get_value('pic_cfl') is not None
        for w in (cb, ui.label_pic_cfl_parcel_fraction,
                  ui.label_pic_cfl_control,
                  ui.lineedit_keyword_pic_cfl_parcel_fraction):
            self.set_widget_enabled(w, enable)
        enable = self.project.get_value('pic_collision_damping', default=False)
        for w in (ui.label_pic_cd_e,
                  ui.lineedit_keyword_pic_cd_e):
            self.set_widget_enabled(w, enable)
        if not enable:
            ui.lineedit_keyword_pic_cd_e.setText('')

    def handle_pic_collision_damping(self):
        ui = self.ui.solids
        enable = self.project.get_value('pic_collision_damping', default=False)
        for w in (ui.label_pic_cd_e,
                  ui.lineedit_keyword_pic_cd_e):
            self.set_widget_enabled(w, enable)
        if not enable:
            self.retain_keyword('pic_cd_e')
            self.unset_keyword('pic_cd_e')
            ui.lineedit_keyword_pic_cd_e.setText('')
        else:
            val = self.get_retained_keyword('pic_cd_e', default=0.8)
            self.update_keyword('pic_cd_e', val)



    def set_pic_defaults(self):
        for (key, default, *rest) in PIC_defaults:
            self.set_keyword_default(key, default)

        #-  Selection enables ‘Solids’ task pane menu
        #-  Selection enables ‘Particle-in-Cell’ task pane menu
        #-  Sets keyword DES_INTERP_ON=.TRUE.
        #-  Sets keyword DES_INTERP_MEAN_FIELDS=.TRUE.
        #-  Sets keyword DES_ONEWAY_COUPLED=.FALSE.
        #-  Sets keyword DES_INTERP_SCHEME=’LINEAR_HAT’
        #-  Sets keyword GENER_PART_CONFIG = .TRUE.

        self.update_keyword('des_interp_on', True)
        self.update_keyword('des_interp_mean_fields', True)
        #self.update_keyword('des_oneway_coupled', False)
        self.unset_keyword('des_oneway_coupled') # False is default
        self.update_keyword('des_interp_scheme', 'LINEAR_HAT')
        self.update_keyword('gener_part_config', True)
        self.unset_keyword('particles') # particles and gener_part_config cannot both be set

    def clear_pic_defaults(self):
        for (key, default, *rest) in PIC_defaults:
            if key != 'des_epg_clip': # EP_STAR/DES_EPG_CLIP is not PIC specific
                self.unset_keyword(key)

        #self.unset_keyword('gener_part_config') # will clobber keyword if model == DEM
        #self.unset_keyword('des_interp_on')
        #self.unset_keyword('des_interp_mean_fields')
        #self.unset_keyword('des_oneway_coupled')

        #self.unset_keyword('des_interp_scheme')
        key = 'des_interp_scheme'
        if self.project.get_value(key) == 'LINEAR_HAT':
            # This setting is required for PIC but not available for other solvers
            self.unset_keyword(key)
