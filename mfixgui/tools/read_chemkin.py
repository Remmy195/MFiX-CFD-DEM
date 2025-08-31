#!/usr/bin/env python

import sys, os
import re

def format_rxn_eqn(rxn_eqn):
    reaction_tmp = rxn_eqn.replace("-", "_")

    M_rxn = False # if there are third body species in the reaction
    M_model = "" # models used for third body reaction
    rev_rxn = False # if it is reversible reaction

    # check for third-body reactions
    if('(+M)' in reaction_tmp):
        reaction_tmp = reaction_tmp.replace('(+M)', '')
        M_rxn = True
        M_model = "M_all"
    elif('+M' in reaction_tmp):
        reaction_tmp = reaction_tmp.replace('+M', '')
        M_rxn = True
        M_model = "M_all"
    elif("(+" in reaction_tmp and ")" in reaction_tmp): # a specific species is used as third body
        M_rxn = True
        index_start = reaction_tmp.index("(+")
        index_end = reaction_tmp.index(")")
        M_model = "M_"+reaction_tmp[index_start+2:index_end]
        reaction_tmp = reaction_tmp.replace(reaction_tmp[index_start:index_end+1], '')

    reaction_tmp = reaction_tmp.replace("(", "_")
    reaction_tmp = reaction_tmp.replace(")", " ")
    # check for reversibale reactions
    if('<=>' in reaction_tmp):
       reaction_tmp = reaction_tmp.replace('<=', '--')
       rev_rxn = True
    elif('=>' in reaction_tmp):
       reaction_tmp = reaction_tmp.replace('=', '--')
    elif('=' in reaction_tmp):
        rev_rxn = True
        reaction_tmp = reaction_tmp.replace('=', '-->')

    reactants, products = reaction_tmp.split('-->')

    reactants = re.sub(r'\s+', '', reactants).split('+')
    products = re.sub(r'\s+', '', products).split('+')

    stoichiometry_rec = []
    species_rec = []
    stoichiometry_prod = []
    species_prod = []
    stoi_rec = {}
    stoi_prod = {}

    for reactant in reactants:
        if re.match(r'(\d+\.\d+|\d+)(.*)', reactant):
            stoichiometry, species = re.match(r'(\d+\.\d+|\d+)(.*)', reactant).groups()
        else:
            stoichiometry = 1.0
            species = reactant
        if(species in species_rec):
            stoichiometry_rec[species_rec.index(species)] += float(stoichiometry)
            stoi_rec[species] += float(stoichiometry)
        else:
            stoichiometry_rec.append(float(stoichiometry))
            species_rec.append(species)
            stoi_rec[species] = float(stoichiometry)

    for product in products:
        if re.match(r'(\d+\.\d+|\d+)(.*)', product):
            stoichiometry, species = re.match(r'(\d+\.\d+|\d+)(.*)', product).groups()
        else:
            stoichiometry = 1.0
            species = product

        if(species in species_prod):
            stoichiometry_prod[species_prod.index(species)] += float(stoichiometry)
            stoi_prod[species] += float(stoichiometry)
        else:
            stoichiometry_prod.append(float(stoichiometry))
            species_prod.append(species)
            stoi_prod[species] = float(stoichiometry)

    # write this reaction as a string
    if(stoichiometry_rec[0] == 1.0):
        rxn_eq = species_rec[0]
    else:
        rxn_eq = str(stoichiometry_rec[0]) + "*" + species_rec[0]
    for rr in range(1, len(species_rec)):
        if(stoichiometry_rec[rr] == 1.0):
            rxn_eq += " + "+ species_rec[rr]
        else:
            rxn_eq += " + "+ str(stoichiometry_rec[rr]) + "*" + species_rec[rr]

    if(stoichiometry_prod[0] == 1.0):
        rxn_eq += " --> " + species_prod[0]
    else:
        rxn_eq += " --> " + str(stoichiometry_prod[0]) + "*" + species_prod[0]
    for pp in range(1, len(species_prod)):
        if(stoichiometry_prod[pp] == 1.0):
            rxn_eq += " + " + species_prod[pp]
        else:
            rxn_eq += " + " + str(stoichiometry_prod[pp]) + "*" + species_prod[pp]

    return rxn_eq, M_rxn, M_model, rev_rxn, stoi_rec, stoi_prod


def read_rate_parameters(rxn_param, species_names, reversible_rxn, third_body_rxn, unit_A, unit_Ea):

    arrhenius_rev = []
    rxn_order_species = []
    rxn_order = []
    rxn_order_species_rev = []
    rxn_order_rev = []
    press_arrhen = []
    press_rxn_type = ""
    press_coeff = []
    press_model_name = "Lindemann" # it will be overwritten if other model is given
    lt_coeff = []
    lt_coeff_rev = []
    fit_coeff = []
    fit_model=""
    M_species = []
    M_coeff = []

    species_names_lower = [a.lower() for a in species_names]
    for ll in range(len(rxn_param)):
        rxn_param_line = rxn_param[ll].strip()
        if (not rxn_param_line.startswith('!')) and (not 'end' in rxn_param_line.lower()):
            # In CHEMKIN format, multiple keywords can be written in one line. BTW, We do not allow this in MFiX.
            # The parameters for each keyword are enclosed by a slash-delimited (/) field.
            # Some keyword may not have any parameters. Therefore, there may not be slash for some keywords.

            # Find the positions of all "/"
            slash_index = [pos for pos, char in enumerate(rxn_param_line) if char == '/']
            if(len(slash_index) %2 != 0.0):
                raise ValueError("Missing '/' in line "+ rxn_param_line +". The parameters for each keyword must be enclosed by a slash-delimited (/) field." )

            # auxiliary keywords in Chemkin
            # DUP or DUPLICATE is used to indicated the duplicate reactions with different parameters
            #     In MFiX, they will be written as different reactions for duplicate reactions
            # REV: arrhenius parameters are given for reverse reactions
            # FORD or RORD: reaction orders for forward and reverse reactions
            # LOW or HIGH: extra arrhnius parameters for low or high pressure.
            #        LOW: for fall-off reactions.
            #        HIGH: for chemically activated bimolecular reactions
            # TROE or SRI: models to calculate the ratefor pressure-dependent reactions, with required coefficients
            # PLOG: pressure dependentce through logarithmic interpolation, with required coefficients
            # LT or RLT: for Landau-Teller reactions (forward or reverse), with required coefficients
            # JAN or FIT1: for optional rate fit expressions, with required coefficients
            # UNITS: units for a specific reaction, which are differed from the default units or the one given in the first line.
            # Species Name: for  neutral third bdy efficiency

            # keywords not supported in MFiX.
            # EXCI: energy loss parameter, in units of electron volts.
            #       This usually used for eletron-impact excitation reactions. We do not support it here.
            #
            # TDEP: for reaction rates depending on specific species temperature.
            #       In MFiX, all species in one phase have the same temperature. So, we do not support it here.
            #
            # CHEB, PCHEB, TCHEB: Chebyshev Polynomial Rate Expressions. This is in ANSYS Chemkin. We do not support it for now.
            # COLLEFF: Efficiency of Collision Frequency Expression. This is in ANSYS Chemkin. We do not support it for now.
            #
            # HV or E can be as reactants or products.
            # When HV is used in the reaction, HV, as a keyword, will be given. But it won't be used.
            #
            # MOME, XSMI: Plasma momentum-transfer collision frequency options

            # Keywords do not need parameters: COLLEFF, DUP, DUPLICATE, MOME, XSMI

            # keywords supported and unsupported in MFiX
            keywords_supported = ['dup', 'duplicate', 'rev', 'ford', 'rord', 'low', 'high', 'troe', 'sri', 'plog', 'lt', 'rlt', \
                                  'jan', 'fit1', 'units']
            keywords_unsupported = ['exci', 'tdep', 'cheb', 'pcheb', 'tcheb', 'colleff', 'hv', 'mome', 'xsmi']

            rxn_param_line_split = rxn_param_line.split('/')
            rxn_param_line_tmp = []
            index_param = 0 # index of keywords and parameters in rxn_param_line_split
            while index_param < len(rxn_param_line_split):
                char = rxn_param_line_split[index_param] # This is a word, not a char - cgw
                if(char != ''):
                    char_lower = char.strip().lower()
                    if not ((char_lower in keywords_supported) or (char_lower in keywords_unsupported) or (char_lower in species_names_lower)):
                        raise ValueError("Unknown keyword "+char+" in line "+rxn_param_line+".")
                    elif(char_lower in keywords_unsupported):
                        raise Warning("Keyword "+char+" in line "+rxn_param_line+" is not supported and will not be used in MFiX.")
                    else:
                        if('dup' in char_lower):
                            index_param += 1 # It is a duplicated reaction. We do not do anything about this. It will be written as another reaction.
                        else:
                            # write all the supported keywords and their parameters into a temporary list
                            rxn_param_line_tmp.append(char_lower)
                            rxn_param_line_tmp.append(rxn_param_line_split[index_param+1])
                            index_param += 2
                else:
                    index_param += 1

            for pos, char in enumerate(rxn_param_line_tmp[::2]):
                parameter_tmp = rxn_param_line_tmp[pos*2+1].split()
                if(char == "rev"):
                    if(reversible_rxn):
                        arrhenius_rev = [float(a) for a in parameter_tmp]
                        if(len(arrhenius_rev)!=3):
                            raise ValueError("3 parameters are required for keyword "+char+", but "+str(len(arrhenius_rev))+" are given in line "+rxn_param_line+".")
                    else:
                        raise ValueError("Arrhenius parameters for the reverse reaction is given for non-reversible reaction "\
                                         "in line "+rxn_param_line+".")
                elif(char == "ford"):
                    for nn in parameter_tmp[::2]:
                        name_tmp = nn.replace("(", "_")
                        name_tmp = name_tmp.replace(")", "")
                        name_tmp = name_tmp.replace("-", "_")
                        rxn_order_species.append(name_tmp)
                    rxn_order += parameter_tmp[1::2]
                elif(char == 'rord'):
                    if(reversible_rxn):
                        for nn in parameter_tmp[::2]:
                            name_tmp = nn.replace("(", "_")
                            name_tmp = name_tmp.replace(")", "")
                            name_tmp = name_tmp.replace("-", "_")
                            rxn_order_species_rev.append(name_tmp)
                        rxn_order_rev += parameter_tmp[1::2]
                    else:
                        raise ValueError("Reaction orders for the reverse reaction is given for non-reversible reaction "\
                                         "in line "+rxn_param_line+".")
                elif(char=="low"):
                    press_arrhen = [float(a) for a in parameter_tmp]
                    press_rxn_type = "falloff"
                    if(len(press_arrhen)!=3):
                        raise ValueError("3 parameters are required for keyword "+char+", but "+str(len(press_arrhen))+ \
                                         " are given in line "+rxn_param_line+".")
                elif(char=="high"):
                    press_arrhen = [float(a) for a in parameter_tmp]
                    press_rxn_type = "bimo"
                    if(len(press_arrhen)!=3):
                        raise ValueError("3 parameters are required for keyword "+char+", but "+str(len(press_arrhen))+ \
                                         " are given in line "+rxn_param_line+".")
                elif(char=='troe'):
                    press_coeff = [float(a) for a in parameter_tmp] # 3 or 4 parameters
                    press_model_name = "Troe"
                    if(len(press_coeff)!=3 and len(press_coeff)!=4):
                        raise ValueError("3 or 4 parameters are required for keyword "+char+", but "+str(len(press_coeff))+ \
                                         " are given in line "+rxn_param_line+".")
                elif(char=='sri'):
                    press_coeff = [float(a) for a in parameter_tmp] # 3 or 5 parameters
                    press_model_name = "SRI"
                    if(len(press_coeff)!=3 and len(press_coeff)!=5):
                        raise ValueError("3 or 5 parameters are required for keyword "+char+", but "+str(len(press_coeff))+ \
                                         " are given in line "+rxn_param_line+".")
                elif(char=='plog'):
                    press_rxn_type = "PLOG"
                    if(len(parameter_tmp)!=4):
                        raise ValueError("4 parameters are required for keyword "+char+", but "+str(len(parameter_tmp))+ \
                                         " are given in line "+rxn_param_line+".")
                    press_arrhen.append([float(a) for a in parameter_tmp[1:]])
                    press_coeff.append(float(parameter_tmp[0]))
                    if(len(press_coeff)>1 and press_coeff[-1] < press_coeff[-2]):
                        raise ValueError("Coefficients in PLOG model must be given in a ascending order of pressure in line "\
			+rxn_param_line+".")
                elif(char=='lt'):
                    lt_coeff = parameter_tmp # 2 parameters
                    if(len(lt_coeff)!=2):
                        raise ValueError("2 parameters are required for keyword "+char+", but "+str(len(lt_coeff))+ \
                                         " are given in line "+rxn_param_line+".")
                elif(char=='rlt'):
                    if(reversible_rxn):
                        lt_coeff_rev = parameter_tmp
                        if(len(lt_coeff_rev)!=2):
                            raise ValueError("2 parameters are required for keyword "+char+", but "+str(len(lt_coeff_rev))+ \
                                             " are given in line "+rxn_param_line+".")
                    else:
                        raise ValueError("LT parameters for the reverse reaction (RTL) is given for non-reversible reaction " \
                                         "in line "+rxn_param_line+".")
                elif(char=='jan'):
                    fit_coeff = parameter_tmp # 9 parameters
                    fit_model = "JAN"
                    if(len(fit_coeff)!=9):
                        raise ValueError("9 parameters are required for keyword "+char+", but "+str(len(fit_coeff))+ \
                                         " are given in line "+rxn_param_line+".")
                elif(char=='fit1'):
                    fit_coeff = parameter_tmp # 4 parameters
                    fit_model = "FIT1"
                    if(len(fit_coeff)!=4):
                        raise ValueError("4 parameters are required for keyword "+char+", but "+str(len(fit_coeff))+ \
                                         " are given in line "+rxn_param_line+".")
                elif(char=='units'): # MOLE(CULE), CAL, KCAL, JOUL, KJOU, KELV(IN), EVOL(TS)
                    for uu in parameter_tmp:
                        uu = uu.lower()
                        if(uu == "cal"):
                            unit_Ea = 4.184
                        elif(uu == "kcal"):
                            unit_Ea = 4184
                        elif(uu == "joul"):
                            unit_Ea = 1
                        elif(uu == "kjou"):
                            unit_Ea = 1000
                        elif(uu == "kelv" or uu == "kelvin"):
                            unit_Ea = 8.3145
                        elif(uu == "evol" or uu == "evolts"):
                            unit_Ea = 96483.04
                        elif(uu == "mole" or uu == "molecule"):
                            unit_A = 1/6.022e23
                        else:
                            raise ValueError("Unsupported units " +uu+" are given in line "+rxn_param_line+".")
                else: # third-body coefficients for species
                    M_species.append(char.upper())
                    if(len(parameter_tmp) >1):
                        raise ValueError("More than one parameters are given for the third-body species, "+char.upper()+ ", in line "+rxn_param_line+".")
                    else:
                        M_coeff.append(parameter_tmp[0])

    if(arrhenius_rev != [] and lt_coeff != [] and lt_coeff_rev==[]):
        raise ValueError("Reverse Landau-Teller parameters must be defined by RLT for reversible reactions with REV.")

    return  (arrhenius_rev, rxn_order_species, rxn_order, rxn_order_species_rev, rxn_order_rev,
             press_rxn_type, press_model_name, press_arrhen, press_coeff,
             lt_coeff, lt_coeff_rev, fit_coeff, fit_model,
             M_species, M_coeff, unit_A, unit_Ea)


def break_string(s, w):
    if len(s) <= w:
        return (s, '')
    else:
        while w > 0 and s[w] != ' ':
            w -= 1
        if w == 0: # no space found, can't break
            return  (s, '')
        s, rest = s[:w], s[w:]
        if len(rest) <= 6: # Prefer a slightly long line to a 'dangler'
            return (s+rest, '')
        else:
            return (s, rest)

def write_string_multiline(line_string):

    string_multiline = ""

    line_string, rest = break_string(line_string, 60)
    if "\"" in line_string:
        quote_sym = "\""
    else:
        quote_sym = ""
    if not rest:
        string_multiline += line_string + " \n"
    else:
        string_multiline += line_string + quote_sym + " & \n"
        while rest:
            line_string, rest = break_string(rest, 60)
            if rest:
                string_multiline += "        " + quote_sym + line_string + "\" & \n"
            else:
                string_multiline += "        " + quote_sym + line_string + " \n"
    return string_multiline


def write_rxn_block(rxn_idx_str, rxn_equ, third_body, reverse_rxn, arrhen, arrhen_rev, \
                    rxn_order_species, rxn_order, rxn_order_species_rev, rxn_order_rev, \
                    press_rxn_type, press_model_name, press_arrhen, press_coeff, \
                    lt_coeff, lt_coeff_rev, fit_coeff, fit_model, \
                    M_species, M_coeff, M_model):
     # write reactions into a string
     rxn_str = "Reaction_"+rxn_idx_str+" {\n"

     eqn_string = write_string_multiline("\""+rxn_equ+"\"")
     rxn_str += "  chem_eq = " + eqn_string

     # add arrhenius coefficients (three numbers, A, beta, E): it MUST be given!
     # for reverse reaction, there are two ways to calculate the reaction rates
     #     reverse_calc = 1: it will be calculated based on the forward reaction constant and equilibrium constant.
     #     reverse_calc = 0: it will be calculated directly based on arrhenius coeffiences for reverse reaction.
     if(reverse_rxn):
         if (arrhen_rev == []):
             rxn_str += "  arrhenius_coeff = " + ' '.join(str("{:.5E}".format(i)) for i in arrhen) +  "\n"
             rxn_str += "  reverse_calc = fromForwardRateConstant \n"
         else:
             rxn_str += "  arrhenius_coeff = " + ' '.join(str("{:.5E}".format(i)) for i in arrhen_rev) +  "\n"
             rxn_str += "  reverse_calc = fromArrheniusCoeff \n"
     # for forward reactions
     else:
         rxn_str += "  arrhenius_coeff = " + ' '.join(str("{:.5E}".format(i)) for i in arrhen) +  "\n"

     # add reaction orders
     if((not reverse_rxn) and (rxn_order_species != [])):
         rxn_order_str = "  rxn_order ="
         for ii, nn in enumerate(rxn_order_species):
             rxn_order_str += " " + nn + ":"+str(rxn_order[ii])
         order_string = write_string_multiline(rxn_order_str)
         rxn_str += order_string
     elif(reverse_rxn and (rxn_order_species_rev != [])):
         rxn_order_str = "  rxn_order ="
         for ii, nn in enumerate(rxn_order_species_rev):
             rxn_order_str += " " + nn + ":"+str(rxn_order_rev[ii])
         order_string = write_string_multiline(rxn_order_str)
         rxn_str += order_string

     # add pressure-related parameters
     if(press_rxn_type != ""):
         if(press_rxn_type == "PLOG"):
             press_str = "  press_rxn_param = " + press_rxn_type
             press_str += " \"press_coeff: " + ' '.join(str("{:.5E}".format(i)) for i in press_coeff)
             press_str +=" arrhenius_press: "
             for ii, pp in enumerate(press_arrhen):
                 press_str += ' '.join(str("{:.5E}".format(i)) for i in pp)
                 if(ii<len(press_arrhen)-1):
                     press_str += ", "
         else:
             press_str = "  press_rxn_param = " + press_model_name + "_" + press_rxn_type
             press_str +=" \"arrhenius_press: " + ' '.join(str("{:.5E}".format(i)) for i in press_arrhen)
             if(press_model_name != "Lindemann"):
                 press_str += " press_coeff: " + ' '.join(str("{:.5E}".format(i)) for i in press_coeff)

         press_str += "\""
         press_string = write_string_multiline(press_str)
         rxn_str += press_string
     # add parameters for Landau-Teller reactions
     if((not reverse_rxn) and lt_coeff != []):
         rxn_str += "  lt_coeff = " + ' '.join(str(i) for i in lt_coeff)
     elif(reverse_rxn and (lt_coeff_rev != [])):
         rxn_str += "  lt_coeff = " + ' '.join(str(i) for i in lt_coeff_rev)

     # add parameters for optional rate fit expressions
     if(fit_coeff != []):
         fit_str = "  fit_coeff = "+ fit_model + "\" coeff: "+ ' '.join(str(i) for i in fit_coeff)
         fit_string = write_string_multiline(fit_str)
         rxn_str += fit_string

     # add parameters for third-body species
     if(third_body):
         third_body_str = "  third_body_param = " + M_model
         if(M_coeff != []):
             third_body_str += " \""
             for ii, pp in enumerate(M_species):
                 name_pp = pp.replace("(", "_")
                 name_pp = name_pp.replace(")", "")
                 name_pp = name_pp.replace("-", "_")
                 third_body_str += name_pp + ":" + str(M_coeff[ii]) + " "
             third_body_str += "\""
         third_body_string = write_string_multiline(third_body_str)
         rxn_str += third_body_string

   #  rxn_str += "  dh = 0 \n  fracdh(0) = 1.0\n "
     rxn_str += "}\n"
     return rxn_str

def write_mfix(mfix_file, mfix_reactions):
    # NOTE: it will write all reactions together, and won't separate homogeneous and heterogeneous reactions.
     with open(mfix_file, 'w') as f:
         for reaction in mfix_reactions:
             f.write(reaction + '\n')





def chemkin_to_mfix(chemkin_lines):
    reactions_lines = []
    index_reactions = []
    rxn_start = False
    species_start = False
    species_list = []
    start_line = [] # List, not line
    for l, line in enumerate(chemkin_lines):
        line_lower = line.strip().lower()
        # Ignore the lines starting with ! (considered as comment)
        if (not line_lower.startswith('!')) and (not 'end' in line_lower):
            # Read species listed in the file
            if('species' in line_lower):
                species_start = True
                rxn_start = False
            # Reaction data must start with the world REACTIONS (or REAC).
            elif ('reac' in line_lower) or ('reactions' in line_lower):
                rxn_start = True
                species_start = False
                start_line = line.split()
            elif species_start:
                species_list += line.split()
            elif rxn_start and ('=' in line):
                # In CHEMKIN, each reaction equation and its Arrhenius coefficients must be contained on one line.
                reactions_lines.append(line.strip())
                index_reactions.append(l)
    # get the unit used in the chemkin mechanism file
    # CAL/MOLE, KCAL/MOLE, JOULES/MOLE, KJOULES/MOLE, KELVINS, or EVOLTS for Ei, and/or MOLES or MOLECULES for Ai
    # T is always in Kelvin.
    # The units written to the mfx file and used in the source code are: CAL/MOLE, MOLES.
    # If no unit is given in the first line of the reaction section, default units (CAL/MOLE, MOLES) will be used.
    unit_fac_Ea = 4.184
    unit_fac_A = 1
    if(len(start_line) > 1):
        start_line_lower = [a.lower() for a in start_line]
        if('cal/mole' in start_line_lower):
            unit_fac_Ea = 4.184
        elif('kcal/mole' in start_line_lower):
            unit_fac_Ea = 4184.0
        elif('joules/mole' in start_line_lower):
            unit_fac_Ea = 1.0
        elif('kjoules/mole' in start_line_lower):
            unit_fac_Ea = 1000.0
        elif('kelvins' in start_line_lower):
            unit_fac_Ea = 8.3145
        elif('evolts' in start_line_lower):
            unit_fac_Ea = 96483.04
        elif('moles' in start_line_lower):
            unit_fac_A = 1
        elif('molecules' in start_line_lower):
            unit_fac_A = 1/6.022e23
        else:
            raise ValueError('Unsupported units "' +' '.join(start_line[1:])+'" found in the first line of reaction block.')

    # getting string for reaction format and reaction-rates related parameters
    mfix_reactions = []

    # number of reactions given in the mechanism
    num_rxn = 0
    for idx, rxn_line in enumerate(reactions_lines):
         num_rxn += 1

	 # split the lines for reaction equations, arrhenius coefficients, and comment (all infos after ! will be considered as comment)
	 # get the index of '!' and the content before '!'
         index_exlc = rxn_line.find('!')
         if(index_exlc !=-1):
             rxn_line = rxn_line[:index_exlc]

         reaction_eq, third_body_rxn, third_body_model, reversible_rxn, stoi_rec, stoi_prod = format_rxn_eqn(''.join(rxn_line.split()[:-3]))
         # arrhenius parameters for this reaction
         arrhenius = [float(a) for a in rxn_line.split()[-3:]]
         # get the index of the next reaction or the end of the file.
         # All data between the line for this reaction and next_index should be the auxiliary information for this reaction.
         if idx <len(reactions_lines)-1:
             next_index = index_reactions[idx+1]
         else:
             next_index = len(chemkin_lines)
         arrhenius_rev, rxn_order_species, rxn_order, rxn_order_species_rev, rxn_order_rev, \
         press_rxn_type, press_model_name, press_arrhen, press_coeff, \
         lt_coeff, lt_coeff_rev, fit_coeff, fit_model, \
         M_species, M_coeff, unit_fac_A_rxn, unit_fac_Ea_rxn \
             = read_rate_parameters(chemkin_lines[index_reactions[idx]+1: next_index], species_list, reversible_rxn, third_body_rxn, unit_fac_A, unit_fac_Ea)

         # get the order of reaction, which is used to convert the units of A from cm-mol to m-kmol
         rxn_order_sum = 0
         rxn_order_sum_rev  = 0
         for ss in list(stoi_rec.keys()):
             if ss in rxn_order_species:
                 rxn_order_sum += float(rxn_order[rxn_order_species.index(ss)])
             else:
                 rxn_order_sum += stoi_rec[ss]
         # modify the factor for Ea to convert from per mol to per kmol.
         unit_fac_Ea_rxn *= 1000
         for ss in list(stoi_prod.keys()):
             if ss in rxn_order_species_rev:
                 rxn_order_sum_rev += float(rxn_order_rev[rxn_order_species_rev.index(ss)])
             else:
                 rxn_order_sum_rev += stoi_prod[ss]

         if(third_body_rxn):
             unit_fac_A_order = (1e-3)**(rxn_order_sum)
             unit_fac_A_order_rev = (1e-3)**(rxn_order_sum_rev)
         else:
             unit_fac_A_order = (1e-3)**(rxn_order_sum-1)
             unit_fac_A_order_rev = (1e-3)**(rxn_order_sum_rev-1)

         arrhenius[0] *= unit_fac_A_rxn*unit_fac_A_order
         arrhenius[2] *= unit_fac_Ea_rxn

         if(press_rxn_type != ""):
             if(press_rxn_type == "PLOG"):
                 for ii in range(len(press_arrhen)):
                     press_arrhen[ii][0] *= unit_fac_A_rxn*unit_fac_A_order
                     press_arrhen[ii][2] *= unit_fac_Ea_rxn
             elif(press_rxn_type == "falloff"):
	     # for fall-off reactions, third body is not required under high pressure.
	     # press_arrhen gives the coefficients for low pressure.
	     # arrhenius gives the coefficients for high pressure.
                 arrhenius[0] /=1e-3
                 press_arrhen[0] *= unit_fac_A_rxn*unit_fac_A_order
                 press_arrhen[2] *= unit_fac_Ea_rxn
             elif(press_rxn_type == "bimo"):
	     # for bimolecular reactions,
	     # press_arrhen gives the coefficients for high pressure.
	     # arrhenius gives the coefficients for low pressure.
                 press_arrhen[0] *= unit_fac_A_rxn*unit_fac_A_order/1e-3
                 press_arrhen[2] *= unit_fac_Ea_rxn

         if(arrhenius_rev != []):
             arrhenius_rev[0] *= unit_fac_A_rxn*unit_fac_A_order_rev
             arrhenius_rev[2] *= unit_fac_Ea_rxn

         reaction_str = write_rxn_block(str(num_rxn), reaction_eq, third_body_rxn, False, \
                                        arrhenius, arrhenius_rev, \
                                        rxn_order_species, rxn_order, \
                                        rxn_order_species_rev, rxn_order_rev, \
                                        press_rxn_type, press_model_name, press_arrhen, press_coeff, \
                                        lt_coeff, lt_coeff_rev, \
                                        fit_coeff, fit_model, \
                                        M_species, M_coeff, third_body_model)

         mfix_reactions.append(reaction_str)

         if reversible_rxn:
             left = reaction_eq.split('-->')[0]
             right = reaction_eq.split('-->')[1]
             reaction_str_rev = write_rxn_block(str(num_rxn)+"_rev", right + ' --> ' + left, third_body_rxn, True, \
                                                arrhenius, arrhenius_rev, \
                                                rxn_order_species, rxn_order, \
                                                rxn_order_species_rev, rxn_order_rev, \
                                                press_rxn_type, press_model_name, press_arrhen, press_coeff, \
                                                lt_coeff, lt_coeff_rev, \
                                                fit_coeff, fit_model, \
                                                M_species, M_coeff, third_body_model)

             mfix_reactions.append(reaction_str_rev)
    return mfix_reactions


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s chemkin_file mfix_file\n"%sys.argv[0])
        sys.exit(-1)

    chemkin_file = sys.argv[1]
    if not os.path.exists(chemkin_file):
        sys.stderr.write("File not found: %s\n" % chemkin_file)
        sys.exit(-1)
    mfix_file =  sys.argv[2]
    if os.path.exists(mfix_file):
        sys.stderr.write("File exists, will not clobber: %s\n" % mfix_file)
        sys.exit(-1)

    with open(chemkin_file, 'r') as f:
        chemkin_lines = f.readlines()
    mfix_reactions = chemkin_to_mfix(chemkin_lines)
    write_mfix(mfix_file, mfix_reactions)



if __name__ == "__main__":
    main()
