import cantera as ct
import sys
import contextlib
import os
import numpy as np


# This is used to suppress the warning about the discontinuity of Cp/R or h/RT.
@contextlib.contextmanager
def suppress_stderr():
    with open(os.devnull, 'w') as fnull:
        stderr_orig = sys.stderr
        sys.stderr = fnull
        try:
            yield
        finally:
            sys.stderr = stderr_orig

def create_cantera_polynomials(species_data, k=True, mu=True):
    # Returns a list, [(keyword, args, value)]
    species_list = []
    ret = []
    if not species_data:
        return ret
    for s, (s_name, s_data) in enumerate(species_data.items(),1):
        temps = s_data['temps']
        if len(temps) > 3:
            raise ValueError("More than two temperature ranges are given for species "
                             +s_name+", which are not supported in Cantera.")
        species_dict = {
            'name': s_name,
            'composition': s_data['composition'],  # Leave empty or infer if needed
            'thermo': {
                'model': 'NASA7',
                'temperature-ranges': [temps[0], temps[1], temps[-1]],
                'data': [s_data['coeffs'][:7], s_data['coeffs'][7:]]
            }
        }
        # TODO: add check for multiple ranges: Cantera does not support multiple ranges.
        # getting the species transport data
        config_transport = s_data.get("config_transport")
        if config_transport in (0,1,2):
            geometry = ("atom", "linear", "nonlinear")[config_transport]
        else:
            raise ValueError("Molecule shape must be 0, 1 or 2, not %s"%config_transport)

        # For Burcat data, we may not have composition information.
        # So, we give the composition based on the molecular weight here.
        if not species_dict['composition']:
            if geometry == "atom":
                species_dict['composition'] = {'H':1} #?
            elif geometry == "linear":
                species_dict['composition'] = {'H':s_data['mol_weight']/1.00794}
            elif geometry == "nonlinear":
                species_dict['composition'] = {'H':s_data['mol_weight']/1.00794}

        species_dict['transport'] = {
            'model': 'gas',
            'geometry': geometry,
            'diameter': s_data["ljsig"],
            'well-depth': s_data["ljeps"],
            'dipole': s_data.get("mu_transport", 0.0),
            'polarizability': s_data.get("alpha_transport", 0.0),
            'rotational-relaxation': s_data.get("zrot_transport", 0.0)
        }
        species = ct.Species.from_dict(species_dict)
        species_list.append(species)

    # cantera mixture based on species information
    with suppress_stderr():
        MIX = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                      transport_model='mixture-averaged',
                      species=species_list,
                      reactions=[])

    for s, s_name in enumerate(species_data.keys(), 1):
        index_mix = MIX.species_index(s_name)
        if k:
            poly_kg = MIX.get_thermal_conductivity_polynomial(index_mix)
            for i in range(len(poly_kg)):
                # Cantera returns np.float64, convert to Python float
                ret.append(("poly_kg", [s,i+1], float(poly_kg[i])))
        if mu:
            poly_mu_g = MIX.get_viscosity_polynomial(index_mix)
            for i in range(len(poly_mu_g)):
                ret.append(("poly_mu_g", [s,i+1], float(poly_mu_g[i])))
    return ret
