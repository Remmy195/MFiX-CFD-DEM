import numpy as np
from qtpy.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QMainWindow,
                            QLabel, QPushButton, QVBoxLayout, QWidget)
from qtpy.QtCore import Qt

from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

from matplotlib.figure import Figure

from mfixgui.tools.cantera_poly import create_cantera_polynomials
from mfixgui.tools.qt import get_icon

def overlaps(r1,r2):
    (a,b),(c,d) = (r1,r2)
    return not (c>b or a>d)
def intersection(r1,r2):
    return [max(r1[0],r2[0]), min(r1[1],r2[1])]

gui = None

plot_titles = ['Specific heat', 'Thermal conductivity', 'Viscosity']
plot_legends = ['cₚ [J/kg.K]', 'k [W/m·K]', 'μ [Pa·s]']

def calc_mu_g(T, MW, ljeps, ljsig):
    # Lennard-Jones model
    Tstar = T/ljeps
    omega_mu = (1.16145/(Tstar**0.14874)
                + 0.52487/np.exp(0.77320*Tstar)
                + 2.16178/np.exp(2.43787*Tstar))
    mu_g = 2.6696e-6 * np.sqrt(MW*T) / (ljsig**2 * omega_mu)
    return mu_g

class PlotWindow(QMainWindow):
    def __init__(self, parent, plot_type):
        global gui
        super().__init__()
        self.parent = parent
        gui = self.parent.parent
        self.plot_type = plot_type # 0: Cp, 1: K, 2: Mu
        self.main = QWidget()
        self.setCentralWidget(self.main)
        vbox = QVBoxLayout(self.main)
        self.canvas = FigureCanvas(Figure(figsize=(8, 6)))
        nt = NavigationToolbar(self.canvas, self)
        self.navigation_toolbar = nt
        self.addToolBar(nt)
        vbox.addWidget(self.canvas)
        hbox = QHBoxLayout()
        vbox.addLayout(hbox)
        b = self.button_prev = QPushButton("Previous species")
        b.setIcon(get_icon('left.svg'))
        b.setFocusPolicy(Qt.NoFocus)
        hbox.addWidget(b)
        b.clicked.connect(self.plot_prev)
        b = self.button_next = QPushButton("Next species")
        b.setLayoutDirection(Qt.RightToLeft)
        b.setIcon(get_icon('right.svg'))
        b.setFocusPolicy(Qt.NoFocus)
        hbox.addWidget(b)
        b.clicked.connect(self.plot_next)
        self.keyReleaseEvent = self.handle_key_event
        self.installEventFilter(self)
        nt.addWidget(QLabel("  "))
        cb = self.checkbox_multiplot = QCheckBox("Multiple plots per page")
        cb.checkable = True
        cb.checked = False
        nt.addWidget(cb)
        cb.clicked.connect(self.toggle_multiplot)
        self.label_multiplot_rows = None
        self.label_multiplot_cols = None
        self.combobox_multiplot_rows = None
        self.combobox_multiplot_cols = None
        self.rows = self.cols = 1
        self.multi = False

    def handle_key_event(self, ev):
        key = ev.key()
        if key in (Qt.Key_Left, Qt.Key_Up, Qt.Key_PageUp):
            self.plot_prev()
        elif key in (Qt.Key_Right, Qt.Key_Down, Qt.Key_PageDown):
            self.plot_next()

    def toggle_multiplot(self, enable):
        nt = self.navigation_toolbar
        if enable:
            self.multi = True
            # Re-adding the same widgets after they have been
            # removed does not work, so we bake 'em fresh every time
            l = self.label_multiplot_rows = QLabel("Rows")
            nt.addWidget(l)
            cb = self.combobox_multiplot_rows = QComboBox()
            cb.activated.connect(self.set_multiplot)
            for i in range(1,11):
                cb.addItem(str(i))
            cb.setCurrentIndex(self.rows-1)
            nt.addWidget(cb)
            l = self.label_multiplot_cols = QLabel("Columns")
            nt.addWidget(l)
            cb = self.combobox_multiplot_cols = QComboBox()
            cb.activated.connect(self.set_multiplot)
            for i in range(1,11):
                cb.addItem(str(i))
            cb.setCurrentIndex(self.cols-1)
            nt.addWidget(cb)
            term = 'page'
        else:
            self.multi = False
            for w in (self.label_multiplot_rows,
                      self.label_multiplot_cols,
                      self.combobox_multiplot_rows,
                      self.combobox_multiplot_cols):
                if w:
                    nt.layout().removeWidget(w)
            term = 'species'

        self.button_prev.setText('Previous ' + term)
        self.button_next.setText('Next ' + term)
        self.do_plot()

    def set_multiplot(self, idx):
        self.rows = self.combobox_multiplot_rows.currentIndex() + 1
        self.cols = self.combobox_multiplot_cols.currentIndex() + 1
        self.do_plot()

    def do_plot(self, species=None):
        if species is None:
            species = self.parent.current_species
        if species is None:
            return
        self.current_species = species
        defined_species = self.parent.defined_species
        n_species = len(defined_species)
        species_list = list(defined_species.keys())
        data = defined_species.get(species)
        t_range = self.parent.t_range

        cvs = self.canvas
        fig = cvs.figure
        fig.clf()

        if self.multi:
            rows = min((n_species+self.cols-1)//self.cols, self.rows) # round up
            cols = min(n_species, self.cols)
        else:
            rows = cols = 1

        axs = fig.subplots(rows, cols, sharex=True)

        get_value = gui.project.get_value
        R = 8.3144598*1000  # Gas constant, J/(K.kmol)

        def plot_cp(ax, species, data):
            temps = data['temps']
            coeffs = data['coeffs']

            MW = data['mol_weight']  # kg/kmol
            r = R/MW  # J/(K.kg)

            def mkplot(ax, a, subrange, lbl):
                T = np.linspace(min(subrange), max(subrange), 500)
                Tpow = [T**i for i in range(5)]
                cp = r * np.dot(a[:5], Tpow)
                ax.plot(T, cp, label=lbl)
                return cp

            if len(temps) == 3:
                tmin, tcom, tmax = temps
                if tmax < tcom:
                    temps = [tmin, tmax]
                    coeffs = coeffs[:7]
                elif tmin > tcom:
                    temps = [tmin, tmax]
                    coeffs = coeffs[7:]

            labels = ["%gK-%gK" % (temps[i], temps[i+1]) for i in range(len(temps)-1)]

            cps = []
            for i in range(len(temps)-1):
                if overlaps(temps[i:i+2], t_range):
                    cp = mkplot(ax, coeffs[i*7: (i+1)*7], intersection(temps[i:i+2], t_range),
                                labels[i])
                    cps.append(cp)

            if t_range[0] < temps[0]:
                ax.plot([t_range[0],temps[0]], [cps[0][0]]*2, 'k--', alpha=.5)
            if t_range[1] > temps[-1]:
                ax.plot([temps[-1],t_range[1]], [cps[-1][-1]]*2, 'k--', alpha=.5)
            ax.legend()

        def plot_k(ax, species, data):
            T = np.linspace(t_range[0], t_range[1], 500)
            kg_model = get_value("kg_model")
            if kg_model == "LENNARD_JONES":
                # This has a lot of overlap with plot_cp
                temps = data['temps']
                coeffs = data['coeffs']
                ljeps = data['ljeps']
                ljsig = data['ljsig']
                MW = data['mol_weight']  # kg/kmol
                r = R/MW  # J/(K.kg)

                def mkplot(ax, a, subrange, lbl):
                    T = np.linspace(min(subrange), max(subrange), 500)
                    mu_g = calc_mu_g(T, MW, ljeps, ljsig)
                    Tpow = [T**i for i in range(5)]
                    cp = np.dot(a[:5], Tpow)
                    kg = (cp+1.25) * r * mu_g
                    ax.plot(T, kg, label=lbl)
                    return kg

                if len(temps) == 3:
                    tmin, tcom, tmax = temps
                    if tmax < tcom:
                        temps = [tmin, tmax]
                        coeffs = coeffs[:7]
                    elif tmin > tcom:
                        temps = [tmin, tmax]
                        coeffs = coeffs[7:]

                labels = ["%gK-%gK" % (temps[i], temps[i+1]) for i in range(len(temps)-1)]

                kgs = []
                for i in range(len(temps)-1):
                    if overlaps(temps[i:i+2], t_range):
                        kg = mkplot(ax, coeffs[i*7: (i+1)*7], intersection(temps[i:i+2], t_range),
                                    labels[i])
                        kgs.append(kg)

                if t_range[0] < temps[0]:
                    ax.plot([t_range[0],temps[0]], [kgs[0][0]]*2, 'k--', alpha=.5)
                if t_range[1] > temps[-1]:
                    ax.plot([temps[-1],t_range[1]], [kgs[-1][-1]]*2, 'k--', alpha=.5)
                ax.legend()


            elif kg_model == "CANTERA_POLY":
                poly_coeffs = np.zeros(5)
                species_index = list(defined_species.keys()).index(species) + 1
                for (key, args, val) in create_cantera_polynomials(defined_species, k=True, mu=False):
                    if args[0] == species_index:
                        poly_coeffs[args[1]-1] = val

                logT = np.log(T)
                logTpow = [logT**i for i in range(5)]
                poly_val = np.dot(poly_coeffs, logTpow)
                kg = np.sqrt(T) * poly_val
                ax.plot(T, kg)
            else:
                raise ValueError("Invalid thermal conductivity model: %s"%kg_model)

        def plot_mu(ax, species, data):
            T = np.linspace(t_range[0], t_range[1], 500)
            mu_g_model = get_value("mu_g_model")
            if mu_g_model == "LENNARD_JONES":
                ljeps = data['ljeps']
                ljsig = data['ljsig']
                MW = data['mol_weight']
                mu_g = calc_mu_g(T, MW, ljeps, ljsig)

            elif mu_g_model == "CANTERA_POLY":
                poly_coeffs = np.zeros(5)
                species_index = list(defined_species.keys()).index(species) + 1
                for (key, args, val) in create_cantera_polynomials(defined_species, k=False, mu=True):
                    if args[0] == species_index:
                        poly_coeffs[args[1]-1] = val
                logT = np.log(T)
                logTpow = [logT**i for i in range(5)]
                poly_val = np.dot(poly_coeffs, logTpow)
                mu_g = np.sqrt(T) * poly_val**2
            else:
                raise ValueError("Invalid viscosity model: %s"%mu_g_model)
            ax.plot(T, mu_g)

        def do_plot(ax, species, data, plot_type):
            ax.grid()
            try:
                if plot_type == 0: # Cp
                    plot_cp(ax, species, data)
                elif self.plot_type == 1: # K
                    plot_k(ax, species, data)
                elif self.plot_type == 2: # Mu
                    plot_mu(ax, species, data)
            except Exception as e:
                gui.error(str(e), popup=True)
                return


        idx = species_list.index(species)
        title = plot_titles[self.plot_type]
        if self.plot_type in (1,2):
            model = get_value('kg_model' if self.plot_type==1 else 'mu_g_model')
            if model == 'LENNARD_JONES':
                title += ' (Lennard-Jones)'
            elif model == 'CANTERA_POLY':
                title += ' (Cantera polynomial)'

        if rows*cols == 1:
            ax = axs
            do_plot(ax, species, data, self.plot_type)
            fig.suptitle(title)
            ax.set_title(species)
            ax.set_xlabel('Temperature [K]')
            ax.set_ylabel(plot_legends[self.plot_type])
            self.setWindowTitle('%s for %s (species %s of %s)' %
                                (plot_titles[self.plot_type], species, idx+1, n_species))
        else:
            # Back up a bit if we're close to the end
            if n_species - idx < rows*cols:
                idx = max(n_species - (rows*cols), 0)
                species = species_list[idx]

            for (i, ax) in enumerate(axs.flat):
                if idx+i >= n_species:
                    break
                species = species_list[idx+i]
                data = defined_species.get(species)
                do_plot(ax, species, data, self.plot_type)
                ax.set_title(species)

            first = idx+1
            last = min(first + n_species - 1, first + rows*cols - 1)
            title1 = plot_titles[self.plot_type]
            # Pluralize properly
            plural = (title1 + 's' if not title1.endswith('y')
                      else title1[:-1]+'ies')
            self.setWindowTitle('%s for species %s - %s of %s' %
                                (plural, first, last, n_species))
            fig.suptitle(title)
            #https://stackoverflow.com/questions/6963035/how-to-set-common-axes-labels-for-subplots
            fig.supxlabel('Temperature [K]')
            fig.supylabel(plot_legends[self.plot_type])

        self.button_prev.setEnabled(
            n_species > 1 and idx > 0)
        self.button_next.setEnabled(
            n_species > 1 and idx < n_species - (rows*cols))

        if not self.isVisible():
            geo = self.parent.geometry()
            x,y = geo.x(), geo.y()
            self.setGeometry(x+20,y+20,800,600)
            self.show()

        fig.set_constrained_layout(True)
        fig.set_tight_layout(False)
        cvs.draw()


    def plot_next(self):
        species = self.current_species
        if not species:
            return
        species_list = list(self.parent.defined_species.keys())
        n_species = len(species_list)
        if species not in species_list:
            return
        idx = species_list.index(species)
        idx0 = idx
        if self.multi:
            idx += self.rows*self.cols
            if idx >= n_species - (self.rows*self.cols):
                idx = n_species - (self.rows*self.cols)
        else:
            idx += 1
            if idx >= n_species:
                return
        if idx == idx0:
            return
        self.do_plot(species_list[idx])
        self.parent.ui.tablewidget_defined_species.selectRow(idx)

    def plot_prev(self):
        species = self.current_species
        if not species:
            return
        species_list = list(self.parent.defined_species.keys())
        n_species = len(species_list)
        if species not in species_list:
            return
        idx = species_list.index(species)
        idx0 = idx
        if self.multi:
            idx -= self.rows*self.cols
            if idx < 0:
                idx = 0
        else:
            idx -= 1
            if idx < 0:
                return
        if idx == idx0:
            return
        self.do_plot(species_list[idx])
        self.parent.ui.tablewidget_defined_species.selectRow(idx)
