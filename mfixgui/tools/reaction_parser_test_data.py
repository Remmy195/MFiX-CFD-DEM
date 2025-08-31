#!/usr/bin/env python

data = """
Phase_Change { chem_eq = "Solid --> Gas" }
Drying { chem_eq = "Liquid --> Vapor" }
Char_Combustion1 { chem_eq = "2FC1 + O2 --> 2CO"}
Char_Combustion2 { chem_eq = "2FC2 + O2 --> 2CO"}
Char_Combustion3 { chem_eq = '2.2FC1 + O2 --> 2CO + 0.2sSoot'}
CH4_Combustion { chem_eq = "CH4 + 2O2 --> CO2 + 2H2O" }
CH4_Catalytic1 { chem_eq = "CH4 + 2O2 + 0.FC1 --> CO2 + 2H2O" }
CH4_Catalytic2 { chem_eq = "CH4 + 2O2 --> CO2 + 2H2O + 0.FC1" }
CH4_Catalytic3 { chem_eq = "CH4 + 2O2 + 0.FC1 --> CO2 + 2H2O + 0.FC1" }
Char_to_Char { chem_eq = 'FC1 --> FC2' }  ! Hot --> Cold
Ash_to_Ash { chem_eq = 'Ash2 --> FlyAsh' } ! Cold --> Hot
Skippy_RxN { chem_eq = 'None' }
Ozone_Decomp { chem_eq = "O3 --> 1.5O2" }
RX1F { chem_eq = "SiH4 --> SiH2 + H2"}
RX1R { chem_eq = "SiH2 + H2 --> SiH4"}
RX2F { chem_eq = "Si2H6 --> SiH4 + SiH2"}
RX2R { chem_eq = "SiH4 + SiH2 --> Si2H6"}
RX3  { chem_eq = "SiH4 --> Si + 2H2"}
RX4  { chem_eq = "SiH2 --> Si + H2"}
CHAR_COMBUSTION { chem_eq = "2FC1 + O2 --> 2CO"}
CHAR_CO2  { chem_eq = "FC1 + CO2 --> 2CO"}           ! Forward
CHAR_CO2r { chem_eq = "2CO + 0.FC1 --> Soot + CO2"}  ! Reverse
CO_Combustion { chem_eq = "CO + 0.5O2 --> CO2" }
Combustion_s1 { chem_eq = "2FC1 + O2 --> 2CO"}
Combustion_s2 { chem_eq = "2FC2 + O2 --> 2CO"}
Char_CO2_s1 { chem_eq = "FC1 + CO2 --> 2CO"}            ! Forward
Char_CO2_s1r { chem_eq = "2CO + 0.FC1 --> Soot + CO2"}  ! Reverse
Char_CO2_s2 { chem_eq = "FC2 + CO2 --> 2CO"}            ! Forward
Char_CO2_s2r { chem_eq = "2CO + 0.FC2 --> Soot + CO2"}  ! Reverse
CO_Combustion1 { chem_eq = "CO + 0.5O2 --> CO2" }
Char_to_Char1 { chem_eq = 'FC2 --> FC1' }
Ash_to_Ash1{ chem_eq = 'Ash2 --> Ash1' }
Drying1 { chem_eq = 'Moisture --> H2O' }
Pyrolysis {
chem_eq = 'Biomass --> ' &
'0.9639 * CO + 0.8771 * CO2 + 0.3491 * CH4 + ' &
'1.6276 * H2 + 1.4210 * H2O'
DH = 3.585d3     ! (cal/mol-biomass)  150.0 J/g-biomass
fracDH(1) = 1.0  ! assign to the coal phase
}
Charring {
chem_eq = 'Biomass --> 8.3334 * Char'
DH = 3.585d3     ! (cal/mol-biomass)  150.0 J/g-biomass
fracDH(1) = 1.0  ! assign to the coal phase
}
Tarring {
chem_eq = 'Biomass --> Tar'
DH = 3.585d3     ! (cal/mol-biomass)  150.0 J/g-biomass
fracDH(1) = 1.0  ! assign to the coal phase
}
AtoR { chem_eq = "A --> 3R" }

Reaction_1 {
  chem_eq = "H2 --> 2.0*H"
  arrhenius_coeff = 4.5770e+19 -1.400 104400.00
  third_body_param = M_all "H2:2.50 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1_rev {
  chem_eq = " 2.0*H --> H2 "
  arrhenius_coeff = 4.5770e+19 -1.400 104400.00
  reverse_calc = fromForwardRateConstant
  third_body_param = M_all "H2:2.50 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_2 {
  chem_eq = "H + CH3 --> CH4"
  arrhenius_coeff = 1.2700e+16 -0.630 383.00
  press_rxn_param = Troe_falloff "arrhenius_press: 2.48e+33" &
        " -4.760 2440.0 press_coeff: 0.7830 74.00 2941. 6964."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_2_rev {
  chem_eq = " CH4 --> H + CH3 "
  arrhenius_coeff = 1.2700e+16 -0.630 383.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 2.48e+33" &
        " -4.760 2440.0 press_coeff: 0.7830 74.00 2941. 6964."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_3 {
  chem_eq = "H + CH4 --> H2 + CH3"
  arrhenius_coeff = 6.1400e+05 2.500 9587.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_3_rev {
  chem_eq = " H2 + CH3 --> H + CH4 "
  arrhenius_coeff = 6.1400e+05 2.500 9587.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_4 {
  chem_eq = "2.0*CH3 --> C2H6"
  arrhenius_coeff = 2.2770e+15 -0.690 174.90
  press_rxn_param = Troe_falloff "arrhenius_press: 8.05e+31" &
        " -3.750 981.6 press_coeff: 0.000 570.0 1.000e+30 1.000e+30"
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_4_rev {
  chem_eq = " C2H6 --> 2.0*CH3 "
  arrhenius_coeff = 2.2770e+15 -0.690 174.90
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 8.05e+31" &
        " -3.750 981.6 press_coeff: 0.000 570.0 1.000e+30 1.000e+30"
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_5 {
  chem_eq = "H + C2H5 --> C2H6"
  arrhenius_coeff = 5.2100e+17 -0.990 1580.00
  press_rxn_param = Troe_falloff "arrhenius_press: 1.99e+41" &
        " -7.080 6685.0 press_coeff: 0.8420 125.0 2219. 6882."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_5_rev {
  chem_eq = " C2H6 --> H + C2H5 "
  arrhenius_coeff = 5.2100e+17 -0.990 1580.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 1.99e+41" &
        " -7.080 6685.0 press_coeff: 0.8420 125.0 2219. 6882."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_6 {
  chem_eq = "H + C2H6 --> H2 + C2H5"
  arrhenius_coeff = 1.1500e+08 1.900 7530.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_6_rev {
  chem_eq = " H2 + C2H5 --> H + C2H6 "
  arrhenius_coeff = 1.1500e+08 1.900 7530.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_7 {
  chem_eq = "CH3 + C2H6 --> CH4 + C2H5"
  arrhenius_coeff = 5.5500e-04 4.720 3231.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_7_rev {
  chem_eq = " CH4 + C2H5 --> CH3 + C2H6 "
  arrhenius_coeff = 5.5500e-04 4.720 3231.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_8 {
  chem_eq = "H + C2H4 --> C2H5"
  arrhenius_coeff = 9.5690e+08 1.463 1355.00
  press_rxn_param = Troe_falloff "arrhenius_press: 1.42e+39" &
        " -6.642 5769.0 press_coeff: -0.5690 299.0 -9147. 152.4"
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_8_rev {
  chem_eq = " C2H5 --> H + C2H4 "
  arrhenius_coeff = 9.5690e+08 1.463 1355.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 1.42e+39" &
        " -6.642 5769.0 press_coeff: -0.5690 299.0 -9147. 152.4"
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_9 {
  chem_eq = "H + C2H5 --> H2 + C2H4"
  arrhenius_coeff = 2.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_9_rev {
  chem_eq = " H2 + C2H4 --> H + C2H5 "
  arrhenius_coeff = 2.0000e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_10 {
  chem_eq = "2.0*C2H4 --> C2H5 + C2H3"
  arrhenius_coeff = 4.8200e+14 0.000 71530.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_10_rev {
  chem_eq = " C2H5 + C2H3 --> 2.0*C2H4 "
  arrhenius_coeff = 4.8200e+14 0.000 71530.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_11 {
  chem_eq = "CH3 + C2H5 --> CH4 + C2H4"
  arrhenius_coeff = 1.1800e+04 2.450 -2921.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_11_rev {
  chem_eq = " CH4 + C2H4 --> CH3 + C2H5 "
  arrhenius_coeff = 1.1800e+04 2.450 -2921.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_12 {
  chem_eq = "2.0*CH3 --> H + C2H5"
  arrhenius_coeff = 4.7400e+12 0.105 10664.30
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 4.740000e+12 1.050000e-01 1.066430e+04," &
        " 2.570000e+13 -9.600000e-02 1.140610e+04, 3.100000e+14" &
        " -3.620000e-01 1.337250e+04, 2.150000e+10 8.850000e-01" &
        " 1.353250e+04, 1.032000e+02 3.230000e+00 1.123610e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_12_rev {
  chem_eq = " H + C2H5 --> 2.0*CH3 "
  arrhenius_coeff = 4.7400e+12 0.105 10664.30
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 4.740000e+12 1.050000e-01 1.066430e+04," &
        " 2.570000e+13 -9.600000e-02 1.140610e+04, 3.100000e+14" &
        " -3.620000e-01 1.337250e+04, 2.150000e+10 8.850000e-01" &
        " 1.353250e+04, 1.032000e+02 3.230000e+00 1.123610e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_13 {
  chem_eq = "H + C2H3 --> C2H4"
  arrhenius_coeff = 6.0800e+12 0.270 280.00
  press_rxn_param = Troe_falloff "arrhenius_press: 1.40e+30" &
        " -3.860 3320.0 press_coeff: 0.7820 207.5 2663. 6095."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_13_rev {
  chem_eq = " C2H4 --> H + C2H3 "
  arrhenius_coeff = 6.0800e+12 0.270 280.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 1.40e+30" &
        " -3.860 3320.0 press_coeff: 0.7820 207.5 2663. 6095."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_14 {
  chem_eq = "C2H4 --> H2 + C2H2"
  arrhenius_coeff = 2.6100e+16 0.000 67823.00
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_14_rev {
  chem_eq = " H2 + C2H2 --> C2H4 "
  arrhenius_coeff = 2.6100e+16 0.000 67823.00
  reverse_calc = fromForwardRateConstant
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_15 {
  chem_eq = "H + C2H4 --> H2 + C2H3"
  arrhenius_coeff = 5.0700e+07 1.930 12950.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_15_rev {
  chem_eq = " H2 + C2H3 --> H + C2H4 "
  arrhenius_coeff = 5.0700e+07 1.930 12950.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_16 {
  chem_eq = "CH3 + C2H4 --> CH4 + C2H3"
  arrhenius_coeff = 9.7600e+02 2.947 15148.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_16_rev {
  chem_eq = " CH4 + C2H3 --> CH3 + C2H4 "
  arrhenius_coeff = 9.7600e+02 2.947 15148.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_17 {
  chem_eq = "CH3 + C2H4 --> CH4 + C2H3"
  arrhenius_coeff = 8.1300e-05 4.417 8835.80
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_17_rev {
  chem_eq = " CH4 + C2H3 --> CH3 + C2H4 "
  arrhenius_coeff = 8.1300e-05 4.417 8835.80
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_18 {
  chem_eq = "H + C2H2 --> C2H3"
  arrhenius_coeff = 1.7100e+10 1.266 2709.00
  press_rxn_param = Troe_falloff "arrhenius_press: 6.35e+31" &
        " -4.664 3780.0 press_coeff: 0.7880 -1.020e+04 1.000e-30"
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_18_rev {
  chem_eq = " C2H3 --> H + C2H2 "
  arrhenius_coeff = 1.7100e+10 1.266 2709.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 6.35e+31" &
        " -4.664 3780.0 press_coeff: 0.7880 -1.020e+04 1.000e-30"
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_19 {
  chem_eq = "H + C2H3 --> H2 + C2H2"
  arrhenius_coeff = 9.6352e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_19_rev {
  chem_eq = " H2 + C2H2 --> H + C2H3 "
  arrhenius_coeff = 9.6352e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_20 {
  chem_eq = "CH3 + C2H3 --> CH4 + C2H2"
  arrhenius_coeff = 3.9200e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_20_rev {
  chem_eq = " CH4 + C2H2 --> CH3 + C2H3 "
  arrhenius_coeff = 3.9200e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_21 {
  chem_eq = "2.0*C2H3 --> C2H4 + C2H2"
  arrhenius_coeff = 9.6000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_21_rev {
  chem_eq = " C2H4 + C2H2 --> 2.0*C2H3 "
  arrhenius_coeff = 9.6000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_22 {
  chem_eq = "H + C2H --> C2H2"
  arrhenius_coeff = 1.0000e+17 -1.000 0.00
  press_rxn_param = Troe_falloff "arrhenius_press: 3.75e+33" &
        " -4.800 1900.0 press_coeff: 0.6460 132.0 1315. 5566."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_22_rev {
  chem_eq = " C2H2 --> H + C2H "
  arrhenius_coeff = 1.0000e+17 -1.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 3.75e+33" &
        " -4.800 1900.0 press_coeff: 0.6460 132.0 1315. 5566."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_23 {
  chem_eq = "H2 + C2H --> H + C2H2"
  arrhenius_coeff = 4.9000e+05 2.500 560.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_23_rev {
  chem_eq = " H + C2H2 --> H2 + C2H "
  arrhenius_coeff = 4.9000e+05 2.500 560.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_24 {
  chem_eq = "C3H8 --> CH3 + C2H5"
  arrhenius_coeff = 1.2900e+37 -5.840 97380.00
  press_rxn_param = Troe_falloff "arrhenius_press: 5.64e+74" &
        " -15.740 98714.0 press_coeff: 0.3100 50.00 3000. 9000."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_24_rev {
  chem_eq = " CH3 + C2H5 --> C3H8 "
  arrhenius_coeff = 1.2900e+37 -5.840 97380.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 5.64e+74" &
        " -15.740 98714.0 press_coeff: 0.3100 50.00 3000. 9000."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_25 {
  chem_eq = "H + NC3H7 --> C3H8"
  arrhenius_coeff = 1.0000e+14 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_25_rev {
  chem_eq = " C3H8 --> H + NC3H7 "
  arrhenius_coeff = 1.0000e+14 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_26 {
  chem_eq = "H + IC3H7 --> C3H8"
  arrhenius_coeff = 1.0000e+14 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_26_rev {
  chem_eq = " C3H8 --> H + IC3H7 "
  arrhenius_coeff = 1.0000e+14 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_27 {
  chem_eq = "C3H8 + IC3H7 --> C3H8 + NC3H7"
  arrhenius_coeff = 3.0000e+10 0.000 12900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_27_rev {
  chem_eq = " C3H8 + NC3H7 --> C3H8 + IC3H7 "
  arrhenius_coeff = 3.0000e+10 0.000 12900.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_28 {
  chem_eq = "H + C3H8 --> H2 + IC3H7"
  arrhenius_coeff = 6.5000e+05 2.400 4471.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_28_rev {
  chem_eq = " H2 + IC3H7 --> H + C3H8 "
  arrhenius_coeff = 6.5000e+05 2.400 4471.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_29 {
  chem_eq = "H + C3H8 --> H2 + NC3H7"
  arrhenius_coeff = 1.7500e+05 2.690 6450.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_29_rev {
  chem_eq = " H2 + NC3H7 --> H + C3H8 "
  arrhenius_coeff = 1.7500e+05 2.690 6450.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_30 {
  chem_eq = "CH3 + C3H8 --> CH4 + IC3H7"
  arrhenius_coeff = 6.4000e+04 2.170 7520.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_30_rev {
  chem_eq = " CH4 + IC3H7 --> CH3 + C3H8 "
  arrhenius_coeff = 6.4000e+04 2.170 7520.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_31 {
  chem_eq = "C2H3 + C3H8 --> C2H4 + IC3H7"
  arrhenius_coeff = 1.0000e+11 0.000 10400.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_31_rev {
  chem_eq = " C2H4 + IC3H7 --> C2H3 + C3H8 "
  arrhenius_coeff = 1.0000e+11 0.000 10400.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_32 {
  chem_eq = "C2H5 + C3H8 --> C2H6 + IC3H7"
  arrhenius_coeff = 1.0000e+11 0.000 10400.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_32_rev {
  chem_eq = " C2H6 + IC3H7 --> C2H5 + C3H8 "
  arrhenius_coeff = 1.0000e+11 0.000 10400.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_33 {
  chem_eq = "C3H8 + C3H5_A --> IC3H7 + C3H6"
  arrhenius_coeff = 7.9400e+11 0.000 16200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_33_rev {
  chem_eq = " IC3H7 + C3H6 --> C3H8 + C3H5_A "
  arrhenius_coeff = 7.9400e+11 0.000 16200.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_34 {
  chem_eq = "CH3 + C3H8 --> CH4 + NC3H7"
  arrhenius_coeff = 9.0400e-01 3.650 7154.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_34_rev {
  chem_eq = " CH4 + NC3H7 --> CH3 + C3H8 "
  arrhenius_coeff = 9.0400e-01 3.650 7154.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_35 {
  chem_eq = "C2H3 + C3H8 --> C2H4 + NC3H7"
  arrhenius_coeff = 1.0000e+11 0.000 10400.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_35_rev {
  chem_eq = " C2H4 + NC3H7 --> C2H3 + C3H8 "
  arrhenius_coeff = 1.0000e+11 0.000 10400.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_36 {
  chem_eq = "C2H5 + C3H8 --> C2H6 + NC3H7"
  arrhenius_coeff = 1.0000e+11 0.000 10400.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_36_rev {
  chem_eq = " C2H6 + NC3H7 --> C2H5 + C3H8 "
  arrhenius_coeff = 1.0000e+11 0.000 10400.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_37 {
  chem_eq = "C3H8 + C3H5_A --> NC3H7 + C3H6"
  arrhenius_coeff = 7.9400e+11 0.000 20500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_37_rev {
  chem_eq = " NC3H7 + C3H6 --> C3H8 + C3H5_A "
  arrhenius_coeff = 7.9400e+11 0.000 20500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_38 {
  chem_eq = "H + IC3H7 --> CH3 + C2H5"
  arrhenius_coeff = 2.0000e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_38_rev {
  chem_eq = " CH3 + C2H5 --> H + IC3H7 "
  arrhenius_coeff = 2.0000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_39 {
  chem_eq = "CH3 + C2H3 --> C3H6"
  arrhenius_coeff = 2.5000e+13 0.000 0.00
  press_rxn_param = Troe_falloff "arrhenius_press: 4.27e+58" &
        " -11.940 9769.8 press_coeff: 0.1750 1341. 6.000e+04" &
        " 1.014e+04"
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_39_rev {
  chem_eq = " C3H6 --> CH3 + C2H3 "
  arrhenius_coeff = 2.5000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 4.27e+58" &
        " -11.940 9769.8 press_coeff: 0.1750 1341. 6.000e+04" &
        " 1.014e+04"
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_40 {
  chem_eq = "CH3 + C2H3 --> H + C3H5_A"
  arrhenius_coeff = 4.1200e+29 -4.950 8000.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 4.120000e+29 -4.950000e+00 8.000000e+03," &
        " 4.860000e+30 -5.030000e+00 1.130000e+04, 5.300000e+29" &
        " -4.570000e+00 1.440000e+04, 1.320000e+30 -4.540000e+00" &
        " 1.930000e+04, 5.160000e+28 -4.030000e+00 2.380000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_40_rev {
  chem_eq = " H + C3H5_A --> CH3 + C2H3 "
  arrhenius_coeff = 4.1200e+29 -4.950 8000.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 4.120000e+29 -4.950000e+00 8.000000e+03," &
        " 4.860000e+30 -5.030000e+00 1.130000e+04, 5.300000e+29" &
        " -4.570000e+00 1.440000e+04, 1.320000e+30 -4.540000e+00" &
        " 1.930000e+04, 5.160000e+28 -4.030000e+00 2.380000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_41 {
  chem_eq = "CH3 + C2H3 --> H + C3H5_A"
  arrhenius_coeff = 5.7300e+15 -0.770 1195.90
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 5.730000e+15 -7.700000e-01 1.195900e+03," &
        " 2.060000e+13 -7.400000e-02 1.428700e+03, 4.480000e+10" &
        " 6.000000e-01 1.421600e+03, 4.100000e+06 1.710000e+00" &
        " 1.056900e+03, 1.370000e-01 3.910000e+00 -3.535500e+02"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_41_rev {
  chem_eq = " H + C3H5_A --> CH3 + C2H3 "
  arrhenius_coeff = 5.7300e+15 -0.770 1195.90
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 5.730000e+15 -7.700000e-01 1.195900e+03," &
        " 2.060000e+13 -7.400000e-02 1.428700e+03, 4.480000e+10" &
        " 6.000000e-01 1.421600e+03, 4.100000e+06 1.710000e+00" &
        " 1.056900e+03, 1.370000e-01 3.910000e+00 -3.535500e+02"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_42 {
  chem_eq = "C3H6 --> CH3 + C2H3"
  arrhenius_coeff = 1.8800e+78 -18.700 130000.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.880000e+78 -1.870000e+01 1.300000e+05," &
        " 8.730000e+76 -1.790000e+01 1.320000e+05, 5.800000e+75" &
        " -1.720000e+01 1.340000e+05, 8.120000e+71 -1.580000e+01" &
        " 1.360000e+05, 2.150000e+64 -1.340000e+01 1.350000e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_42_rev {
  chem_eq = " CH3 + C2H3 --> C3H6 "
  arrhenius_coeff = 1.8800e+78 -18.700 130000.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.880000e+78 -1.870000e+01 1.300000e+05," &
        " 8.730000e+76 -1.790000e+01 1.320000e+05, 5.800000e+75" &
        " -1.720000e+01 1.340000e+05, 8.120000e+71 -1.580000e+01" &
        " 1.360000e+05, 2.150000e+64 -1.340000e+01 1.350000e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_43 {
  chem_eq = "C3H6 --> CH3 + C2H3"
  arrhenius_coeff = 1.6900e+59 -13.600 113290.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.690000e+59 -1.360000e+01 1.132900e+05," &
        " 2.000000e+60 -1.370000e+01 1.148900e+05, 6.700000e+54" &
        " -1.180000e+01 1.138400e+05, 1.060000e+47 -9.270000e+00" &
        " 1.115100e+05, 7.290000e+38 -6.700000e+00 1.087400e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_43_rev {
  chem_eq = " CH3 + C2H3 --> C3H6 "
  arrhenius_coeff = 1.6900e+59 -13.600 113290.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.690000e+59 -1.360000e+01 1.132900e+05," &
        " 2.000000e+60 -1.370000e+01 1.148900e+05, 6.700000e+54" &
        " -1.180000e+01 1.138400e+05, 1.060000e+47 -9.270000e+00" &
        " 1.115100e+05, 7.290000e+38 -6.700000e+00 1.087400e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_44 {
  chem_eq = "C3H6 --> H + C3H5_A"
  arrhenius_coeff = 9.1600e+74 -17.600 120000.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 9.160000e+74 -1.760000e+01 1.200000e+05," &
        " 1.730000e+70 -1.600000e+01 1.200000e+05, 1.080000e+71" &
        " -1.590000e+01 1.248600e+05, 6.400000e+65 -1.420000e+01" &
        " 1.250000e+05, 8.050000e+56 -1.150000e+01 1.220000e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_44_rev {
  chem_eq = " H + C3H5_A --> C3H6 "
  arrhenius_coeff = 9.1600e+74 -17.600 120000.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 9.160000e+74 -1.760000e+01 1.200000e+05," &
        " 1.730000e+70 -1.600000e+01 1.200000e+05, 1.080000e+71" &
        " -1.590000e+01 1.248600e+05, 6.400000e+65 -1.420000e+01" &
        " 1.250000e+05, 8.050000e+56 -1.150000e+01 1.220000e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_45 {
  chem_eq = "C3H6 --> H + C3H5_A"
  arrhenius_coeff = 2.9800e+54 -12.300 101200.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.980000e+54 -1.230000e+01 1.012000e+05," &
        " 1.370000e+43 -8.870000e+00 9.636500e+04, 6.280000e+42" &
        " -8.510000e+00 9.800400e+04, 4.730000e+35 -6.260000e+00" &
        " 9.564400e+04, 4.340000e+28 -4.060000e+00 9.311400e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_45_rev {
  chem_eq = " H + C3H5_A --> C3H6 "
  arrhenius_coeff = 2.9800e+54 -12.300 101200.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.980000e+54 -1.230000e+01 1.012000e+05," &
        " 1.370000e+43 -8.870000e+00 9.636500e+04, 6.280000e+42" &
        " -8.510000e+00 9.800400e+04, 4.730000e+35 -6.260000e+00" &
        " 9.564400e+04, 4.340000e+28 -4.060000e+00 9.311400e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_46 {
  chem_eq = "H + C3H5_T --> C3H6"
  arrhenius_coeff = 4.9600e+60 -15.200 18000.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 4.960000e+60 -1.520000e+01 1.800000e+04," &
        " 3.200000e+62 -1.510000e+01 2.010000e+04, 2.310000e+60" &
        " -1.400000e+01 2.190000e+04, 3.690000e+54 -1.200000e+01" &
        " 2.210000e+04, 1.150000e+50 -1.040000e+01 2.330000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_46_rev {
  chem_eq = " C3H6 --> H + C3H5_T "
  arrhenius_coeff = 4.9600e+60 -15.200 18000.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 4.960000e+60 -1.520000e+01 1.800000e+04," &
        " 3.200000e+62 -1.510000e+01 2.010000e+04, 2.310000e+60" &
        " -1.400000e+01 2.190000e+04, 3.690000e+54 -1.200000e+01" &
        " 2.210000e+04, 1.150000e+50 -1.040000e+01 2.330000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_47 {
  chem_eq = "H + C3H5_T --> C3H6"
  arrhenius_coeff = 1.4900e+48 -12.000 7203.30
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.490000e+48 -1.200000e+01 7.203300e+03," &
        " 6.760000e+46 -1.110000e+01 7.629900e+03, 1.090000e+40" &
        " -8.660000e+00 6.447800e+03, 2.380000e+31 -5.730000e+00" &
        " 4.506000e+03, 5.690000e+25 -3.830000e+00 3.250400e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_47_rev {
  chem_eq = " C3H6 --> H + C3H5_T "
  arrhenius_coeff = 1.4900e+48 -12.000 7203.30
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.490000e+48 -1.200000e+01 7.203300e+03," &
        " 6.760000e+46 -1.110000e+01 7.629900e+03, 1.090000e+40" &
        " -8.660000e+00 6.447800e+03, 2.380000e+31 -5.730000e+00" &
        " 4.506000e+03, 5.690000e+25 -3.830000e+00 3.250400e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_48 {
  chem_eq = "H + C3H5_T --> H + C3H5_A"
  arrhenius_coeff = 2.1100e+17 -1.080 1290.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.110000e+17 -1.080000e+00 1.290000e+03," &
        " 9.050000e+29 -4.910000e+00 8.540000e+03, 2.980000e+30" &
        " -4.790000e+00 1.200000e+04, 8.220000e+28 -4.140000e+00" &
        " 1.540000e+04, 2.280000e+29 -4.120000e+00 2.090000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_48_rev {
  chem_eq = " H + C3H5_A --> H + C3H5_T "
  arrhenius_coeff = 2.1100e+17 -1.080 1290.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.110000e+17 -1.080000e+00 1.290000e+03," &
        " 9.050000e+29 -4.910000e+00 8.540000e+03, 2.980000e+30" &
        " -4.790000e+00 1.200000e+04, 8.220000e+28 -4.140000e+00" &
        " 1.540000e+04, 2.280000e+29 -4.120000e+00 2.090000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_49 {
  chem_eq = "H + C3H5_T --> H + C3H5_A"
  arrhenius_coeff = 6.4100e+03 2.610 -3778.40
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 6.410000e+03 2.610000e+00 -3.778400e+03," &
        " 5.190000e+14 -3.000000e-01 1.090400e+03, 8.170000e+11" &
        " 4.900000e-01 1.184600e+03, 2.790000e+09 1.090000e+00" &
        " 1.187500e+03, 6.750000e+03 2.700000e+00 3.738000e+02"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_49_rev {
  chem_eq = " H + C3H5_A --> H + C3H5_T "
  arrhenius_coeff = 6.4100e+03 2.610 -3778.40
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 6.410000e+03 2.610000e+00 -3.778400e+03," &
        " 5.190000e+14 -3.000000e-01 1.090400e+03, 8.170000e+11" &
        " 4.900000e-01 1.184600e+03, 2.790000e+09 1.090000e+00" &
        " 1.187500e+03, 6.750000e+03 2.700000e+00 3.738000e+02"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_50 {
  chem_eq = "H + C3H5_T --> CH3 + C2H3"
  arrhenius_coeff = 3.3100e+16 -0.690 5200.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 3.310000e+16 -6.900000e-01 5.200000e+03," &
        " 9.040000e+16 -8.100000e-01 4.800000e+03, 2.010000e+24" &
        " -2.860000e+00 1.090000e+04, 2.750000e+26 -3.310000e+00" &
        " 1.580000e+04, 3.150000e+32 -4.830000e+00 2.600000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_50_rev {
  chem_eq = " CH3 + C2H3 --> H + C3H5_T "
  arrhenius_coeff = 3.3100e+16 -0.690 5200.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 3.310000e+16 -6.900000e-01 5.200000e+03," &
        " 9.040000e+16 -8.100000e-01 4.800000e+03, 2.010000e+24" &
        " -2.860000e+00 1.090000e+04, 2.750000e+26 -3.310000e+00" &
        " 1.580000e+04, 3.150000e+32 -4.830000e+00 2.600000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_51 {
  chem_eq = "H + C3H5_T --> CH3 + C2H3"
  arrhenius_coeff = 8.0400e+13 -0.140 1150.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 8.040000e+13 -1.400000e-01 1.150000e+03," &
        " 7.170000e+10 6.700000e-01 6.738000e+02, 9.970000e+08" &
        " 1.360000e+00 1.596400e+03, 7.410000e+07 1.570000e+00" &
        " 2.108800e+03, 2.700000e+12 3.200000e-01 6.791800e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_51_rev {
  chem_eq = " CH3 + C2H3 --> H + C3H5_T "
  arrhenius_coeff = 8.0400e+13 -0.140 1150.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 8.040000e+13 -1.400000e-01 1.150000e+03," &
        " 7.170000e+10 6.700000e-01 6.738000e+02, 9.970000e+08" &
        " 1.360000e+00 1.596400e+03, 7.410000e+07 1.570000e+00" &
        " 2.108800e+03, 2.700000e+12 3.200000e-01 6.791800e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_52 {
  chem_eq = "C3H6 --> H + C3H5_S"
  arrhenius_coeff = 7.7100e+69 -16.090 140000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_52_rev {
  chem_eq = " H + C3H5_S --> C3H6 "
  arrhenius_coeff = 7.7100e+69 -16.090 140000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_53 {
  chem_eq = "H + C3H6 --> H2 + C3H5_A"
  arrhenius_coeff = 3.6440e+05 2.455 4361.20
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_53_rev {
  chem_eq = " H2 + C3H5_A --> H + C3H6 "
  arrhenius_coeff = 3.6440e+05 2.455 4361.20
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_54 {
  chem_eq = "CH3 + C3H6 --> CH4 + C3H5_A"
  arrhenius_coeff = 2.2100e+00 3.500 5675.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_54_rev {
  chem_eq = " CH4 + C3H5_A --> CH3 + C3H6 "
  arrhenius_coeff = 2.2100e+00 3.500 5675.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_55 {
  chem_eq = "C2H5 + C3H6 --> C2H6 + C3H5_A"
  arrhenius_coeff = 1.0000e+11 0.000 9800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_55_rev {
  chem_eq = " C2H6 + C3H5_A --> C2H5 + C3H6 "
  arrhenius_coeff = 1.0000e+11 0.000 9800.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_56 {
  chem_eq = "H + C3H6 --> H2 + C3H5_T"
  arrhenius_coeff = 1.4980e+02 3.381 8909.50
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_56_rev {
  chem_eq = " H2 + C3H5_T --> H + C3H6 "
  arrhenius_coeff = 1.4980e+02 3.381 8909.50
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_57 {
  chem_eq = "CH3 + C3H6 --> CH4 + C3H5_T"
  arrhenius_coeff = 8.4000e-01 3.500 11660.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_57_rev {
  chem_eq = " CH4 + C3H5_T --> CH3 + C3H6 "
  arrhenius_coeff = 8.4000e-01 3.500 11660.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_58 {
  chem_eq = "H + C3H6 --> H2 + C3H5_S"
  arrhenius_coeff = 5.1010e+02 3.234 12357.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_58_rev {
  chem_eq = " H2 + C3H5_S --> H + C3H6 "
  arrhenius_coeff = 5.1010e+02 3.234 12357.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_59 {
  chem_eq = "H + C3H6 --> H2 + C3H5_S"
  arrhenius_coeff = 3.9690e+02 3.252 12007.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_59_rev {
  chem_eq = " H2 + C3H5_S --> H + C3H6 "
  arrhenius_coeff = 3.9690e+02 3.252 12007.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_60 {
  chem_eq = "CH3 + C3H6 --> CH4 + C3H5_S"
  arrhenius_coeff = 1.3480e+00 3.500 12850.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_60_rev {
  chem_eq = " CH4 + C3H5_S --> CH3 + C3H6 "
  arrhenius_coeff = 1.3480e+00 3.500 12850.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_61 {
  chem_eq = "H + C3H6 --> NC3H7"
  arrhenius_coeff = 1.0000e+00 1.000 0.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 7.990000e+81 -2.316100e+01 2.223900e+04," &
        " 4.240000e+68 -1.842700e+01 1.966500e+04, 1.040000e+49" &
        " -1.150000e+01 1.535900e+04, 6.200000e+41 -8.892000e+00" &
        " 1.463700e+04, 1.000000e-10 0.000000e+00 0.000000e+00"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_61_rev {
  chem_eq = " NC3H7 --> H + C3H6 "
  arrhenius_coeff = 1.0000e+00 1.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 7.990000e+81 -2.316100e+01 2.223900e+04," &
        " 4.240000e+68 -1.842700e+01 1.966500e+04, 1.040000e+49" &
        " -1.150000e+01 1.535900e+04, 6.200000e+41 -8.892000e+00" &
        " 1.463700e+04, 1.000000e-10 0.000000e+00 0.000000e+00"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_62 {
  chem_eq = "H + C3H6 --> NC3H7"
  arrhenius_coeff = 1.0000e+00 1.000 0.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.850000e+26 -5.830000e+00 3.865800e+03," &
        " 2.820000e+30 -6.490000e+00 5.470800e+03, 3.780000e+28" &
        " -5.570000e+00 5.625100e+03, 1.460000e+25 -4.280000e+00" &
        " 5.247800e+03, 4.220000e+27 -4.390000e+00 9.345800e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_62_rev {
  chem_eq = " NC3H7 --> H + C3H6 "
  arrhenius_coeff = 1.0000e+00 1.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.850000e+26 -5.830000e+00 3.865800e+03," &
        " 2.820000e+30 -6.490000e+00 5.470800e+03, 3.780000e+28" &
        " -5.570000e+00 5.625100e+03, 1.460000e+25 -4.280000e+00" &
        " 5.247800e+03, 4.220000e+27 -4.390000e+00 9.345800e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_63 {
  chem_eq = "H + C3H6 --> CH3 + C2H4"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.540000e+09 1.350000e+00 2.542000e+03," &
        " 7.880000e+10 8.700000e-01 3.599600e+03, 2.670000e+12" &
        " 4.700000e-01 5.431100e+03, 9.250000e+22 -2.600000e+00" &
        " 1.289800e+04, 1.320000e+23 -2.420000e+00 1.650000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_63_rev {
  chem_eq = " CH3 + C2H4 --> H + C3H6 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.540000e+09 1.350000e+00 2.542000e+03," &
        " 7.880000e+10 8.700000e-01 3.599600e+03, 2.670000e+12" &
        " 4.700000e-01 5.431100e+03, 9.250000e+22 -2.600000e+00" &
        " 1.289800e+04, 1.320000e+23 -2.420000e+00 1.650000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_64 {
  chem_eq = "H + C3H6 --> CH3 + C2H4"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.000000e-10 0.000000e+00 0.000000e+00," &
        " 1.000000e-10 0.000000e+00 0.000000e+00, 1.000000e-10" &
        " 0.000000e+00 0.000000e+00, 1.240000e+05 2.520000e+00" &
        " 3.679100e+03, 2.510000e+03 2.910000e+00 3.980900e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_64_rev {
  chem_eq = " CH3 + C2H4 --> H + C3H6 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.000000e-10 0.000000e+00 0.000000e+00," &
        " 1.000000e-10 0.000000e+00 0.000000e+00, 1.000000e-10" &
        " 0.000000e+00 0.000000e+00, 1.240000e+05 2.520000e+00" &
        " 3.679100e+03, 2.510000e+03 2.910000e+00 3.980900e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_65 {
  chem_eq = "H + C3H6 --> IC3H7"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.350000e+44 -1.068000e+01 8.196400e+03," &
        " 2.110000e+57 -1.423000e+01 1.514700e+04, 3.260000e+61" &
        " -1.494000e+01 2.016100e+04, 5.300000e+56 -1.312000e+01" &
        " 2.066700e+04, 1.110000e+50 -1.080000e+01 2.020200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_65_rev {
  chem_eq = " IC3H7 --> H + C3H6 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.350000e+44 -1.068000e+01 8.196400e+03," &
        " 2.110000e+57 -1.423000e+01 1.514700e+04, 3.260000e+61" &
        " -1.494000e+01 2.016100e+04, 5.300000e+56 -1.312000e+01" &
        " 2.066700e+04, 1.110000e+50 -1.080000e+01 2.020200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_66 {
  chem_eq = "H + C3H6 --> IC3H7"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.170000e+130 -3.258000e+01 1.361400e+05," &
        " 2.250000e+29 -5.840000e+00 4.241900e+03, 1.060000e+30" &
        " -5.630000e+00 5.613400e+03, 6.110000e+26 -4.440000e+00" &
        " 5.182300e+03, 2.730000e+23 -3.260000e+00 4.597000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_66_rev {
  chem_eq = " IC3H7 --> H + C3H6 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.170000e+130 -3.258000e+01 1.361400e+05," &
        " 2.250000e+29 -5.840000e+00 4.241900e+03, 1.060000e+30" &
        " -5.630000e+00 5.613400e+03, 6.110000e+26 -4.440000e+00" &
        " 5.182300e+03, 2.730000e+23 -3.260000e+00 4.597000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_67 {
  chem_eq = "CH3 + C2H4 --> NC3H7"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 8.670000e+48 -1.254000e+01 1.820600e+04," &
        " 1.060000e+49 -1.204000e+01 2.000100e+04, 7.670000e+47" &
        " -1.117000e+01 2.236600e+04, 1.810000e+45 -1.003000e+01" &
        " 2.376900e+04, 2.040000e+40 -8.250000e+00 2.421400e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_67_rev {
  chem_eq = " NC3H7 --> CH3 + C2H4 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 8.670000e+48 -1.254000e+01 1.820600e+04," &
        " 1.060000e+49 -1.204000e+01 2.000100e+04, 7.670000e+47" &
        " -1.117000e+01 2.236600e+04, 1.810000e+45 -1.003000e+01" &
        " 2.376900e+04, 2.040000e+40 -8.250000e+00 2.421400e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_68 {
  chem_eq = "CH3 + C2H4 --> NC3H7"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.120000e+43 -1.130000e+01 1.308000e+04," &
        " 7.280000e+39 -9.880000e+00 1.316400e+04, 2.600000e+33" &
        " -7.460000e+00 1.241600e+04, 3.850000e+27 -5.380000e+00" &
        " 1.145500e+04, 1.660000e+21 -3.170000e+00 1.024100e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_68_rev {
  chem_eq = " NC3H7 --> CH3 + C2H4 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.120000e+43 -1.130000e+01 1.308000e+04," &
        " 7.280000e+39 -9.880000e+00 1.316400e+04, 2.600000e+33" &
        " -7.460000e+00 1.241600e+04, 3.850000e+27 -5.380000e+00" &
        " 1.145500e+04, 1.660000e+21 -3.170000e+00 1.024100e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_69 {
  chem_eq = "H + C3H5_A --> H2 + C3H4_A"
  arrhenius_coeff = 1.2320e+03 3.035 2582.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_69_rev {
  chem_eq = " H2 + C3H4_A --> H + C3H5_A "
  arrhenius_coeff = 1.2320e+03 3.035 2582.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_70 {
  chem_eq = "CH3 + C3H5_A --> CH4 + C3H4_A"
  arrhenius_coeff = 3.0000e+12 -0.320 -131.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_70_rev {
  chem_eq = " CH4 + C3H4_A --> CH3 + C3H5_A "
  arrhenius_coeff = 3.0000e+12 -0.320 -131.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_71 {
  chem_eq = "C2H5 + C3H5_A --> C2H6 + C3H4_A"
  arrhenius_coeff = 4.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_71_rev {
  chem_eq = " C2H6 + C3H4_A --> C2H5 + C3H5_A "
  arrhenius_coeff = 4.0000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_72 {
  chem_eq = "C2H3 + C3H5_A --> C2H4 + C3H4_A"
  arrhenius_coeff = 1.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_72_rev {
  chem_eq = " C2H4 + C3H4_A --> C2H3 + C3H5_A "
  arrhenius_coeff = 1.0000e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_73 {
  chem_eq = "C2H2 + C3H5_A --> C2H4 + C3H3"
  arrhenius_coeff = 9.4500e+08 1.430 29879.55
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 9.450000e+08" &
        " 1.430000e+00 2.987955e+04, 3.010000e+09 1.300000e+00" &
        " 3.029287e+04, 8.630000e+17 -1.010000e+00 3.800488e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_74 {
  chem_eq = "2.0*C3H4_A --> C3H5_A + C3H3"
  arrhenius_coeff = 5.0000e+14 0.000 64746.70
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_74_rev {
  chem_eq = " C3H5_A + C3H3 --> 2.0*C3H4_A "
  arrhenius_coeff = 5.0000e+14 0.000 64746.70
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_75 {
  chem_eq = "H + C3H5_S --> H2 + C3H4_A"
  arrhenius_coeff = 3.3330e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_75_rev {
  chem_eq = " H2 + C3H4_A --> H + C3H5_S "
  arrhenius_coeff = 3.3330e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_76 {
  chem_eq = "CH3 + C3H5_S --> CH4 + C3H4_A"
  arrhenius_coeff = 1.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_76_rev {
  chem_eq = " CH4 + C3H4_A --> CH3 + C3H5_S "
  arrhenius_coeff = 1.0000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_77 {
  chem_eq = "H + C3H5_S --> H2 + C3H4_P"
  arrhenius_coeff = 3.3400e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_77_rev {
  chem_eq = " H2 + C3H4_P --> H + C3H5_S "
  arrhenius_coeff = 3.3400e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_78 {
  chem_eq = "CH3 + C3H5_S --> CH4 + C3H4_P"
  arrhenius_coeff = 1.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_78_rev {
  chem_eq = " CH4 + C3H4_P --> CH3 + C3H5_S "
  arrhenius_coeff = 1.0000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_79 {
  chem_eq = "H + C3H5_T --> H2 + C3H4_P"
  arrhenius_coeff = 3.3400e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_79_rev {
  chem_eq = " H2 + C3H4_P --> H + C3H5_T "
  arrhenius_coeff = 3.3400e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_80 {
  chem_eq = "CH3 + C3H5_T --> CH4 + C3H4_P"
  arrhenius_coeff = 1.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_80_rev {
  chem_eq = " CH4 + C3H4_P --> CH3 + C3H5_T "
  arrhenius_coeff = 1.0000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_81 {
  chem_eq = "2.0*C3H5_A --> C3H6 + C3H4_A"
  arrhenius_coeff = 9.5500e+40 -9.300 12470.00
  press_rxn_param = PLOG "press_coeff: 1.000000e+00" &
        " 4.000000e+00 1.000000e+01 arrhenius_press: 4.770000e+40" &
        " -9.300000e+00 1.247000e+04, 3.970000e+32 -6.800000e+00" &
        " 9.180000e+03, 1.460000e+28 -5.500000e+00 7.410000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_81_rev {
  chem_eq = " C3H6 + C3H4_A --> 2.0*C3H5_A "
  arrhenius_coeff = 9.5500e+40 -9.300 12470.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e+00" &
        " 4.000000e+00 1.000000e+01 arrhenius_press: 4.770000e+40" &
        " -9.300000e+00 1.247000e+04, 3.970000e+32 -6.800000e+00" &
        " 9.180000e+03, 1.460000e+28 -5.500000e+00 7.410000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_82 {
  chem_eq = "C2H5 + C3H5_A --> C2H4 + C3H6"
  arrhenius_coeff = 4.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_82_rev {
  chem_eq = " C2H4 + C3H6 --> C2H5 + C3H5_A "
  arrhenius_coeff = 4.0000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_83 {
  chem_eq = "C3H4_P + C3H3 --> C3H4_A + C3H3"
  arrhenius_coeff = 6.1400e+06 1.740 10450.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_83_rev {
  chem_eq = " C3H4_A + C3H3 --> C3H4_P + C3H3 "
  arrhenius_coeff = 6.1400e+06 1.740 10450.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_84 {
  chem_eq = "H + C3H4_P --> H2 + C3H3"
  arrhenius_coeff = 3.5720e+04 2.825 4821.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_84_rev {
  chem_eq = " H2 + C3H3 --> H + C3H4_P "
  arrhenius_coeff = 3.5720e+04 2.825 4821.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_85 {
  chem_eq = "CH3 + C3H4_P --> CH4 + C3H3"
  arrhenius_coeff = 1.8000e+12 0.000 7700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_85_rev {
  chem_eq = " CH4 + C3H3 --> CH3 + C3H4_P "
  arrhenius_coeff = 1.8000e+12 0.000 7700.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_86 {
  chem_eq = "C2H + C3H4_P --> C2H2 + C3H3"
  arrhenius_coeff = 1.0000e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_86_rev {
  chem_eq = " C2H2 + C3H3 --> C2H + C3H4_P "
  arrhenius_coeff = 1.0000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_87 {
  chem_eq = "C2H3 + C3H4_P --> C2H4 + C3H3"
  arrhenius_coeff = 3.3100e+15 -1.440 17273.69
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 3.310000e+15" &
        " -1.440000e+00 1.727369e+04, 1.430000e+20 -2.520000e+00" &
        " 2.747363e+04, 1.000000e+16 -1.010000e+00 4.136310e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_88 {
  chem_eq = "C2H3 + C3H4_P --> CH3 + C4H4"
  arrhenius_coeff = 3.3400e+08 1.250 7662.70
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 3.340000e+08" &
        " 1.250000e+00 7.662700e+03, 6.280000e+16 -1.010000e+00" &
        " 1.504912e+04, 4.440000e+16 -7.800000e-01 2.102974e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_89 {
  chem_eq = "C3H5_A + C3H4_P --> C3H6 + C3H3"
  arrhenius_coeff = 3.0000e+12 0.000 7700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_89_rev {
  chem_eq = " C3H6 + C3H3 --> C3H5_A + C3H4_P "
  arrhenius_coeff = 3.0000e+12 0.000 7700.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_90 {
  chem_eq = "H + C3H4_A --> H2 + C3H3"
  arrhenius_coeff = 6.6250e+03 3.095 5522.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_90_rev {
  chem_eq = " H2 + C3H3 --> H + C3H4_A "
  arrhenius_coeff = 6.6250e+03 3.095 5522.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_91 {
  chem_eq = "CH3 + C3H4_A --> CH4 + C3H3"
  arrhenius_coeff = 1.3000e+12 0.000 7700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_91_rev {
  chem_eq = " CH4 + C3H3 --> CH3 + C3H4_A "
  arrhenius_coeff = 1.3000e+12 0.000 7700.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_92 {
  chem_eq = "C3H5_A + C3H4_A --> C3H6 + C3H3"
  arrhenius_coeff = 2.0000e+11 0.000 7700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_92_rev {
  chem_eq = " C3H6 + C3H3 --> C3H5_A + C3H4_A "
  arrhenius_coeff = 2.0000e+11 0.000 7700.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_93 {
  chem_eq = "C2H3 + C3H4_A --> C2H4 + C3H3"
  arrhenius_coeff = 1.0300e+27 -4.130 17226.03
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 1.030000e+27" &
        " -4.130000e+00 1.722603e+04, 6.330000e+26 -3.860000e+00" &
        " 2.420763e+04, 2.480000e+14 -1.000000e-01 3.006240e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_94 {
  chem_eq = "C2H3 + C3H4_A --> CH3 + C4H4"
  arrhenius_coeff = 1.2500e+11 0.310 17562.92
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 1.250000e+11" &
        " 3.100000e-01 1.756292e+04, 1.700000e+17 -1.220000e+00" &
        " 2.749888e+04, 2.930000e+01 3.470000e+00 2.912062e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_95 {
  chem_eq = "H + C3H4_A --> H + C3H4_P"
  arrhenius_coeff = 2.4400e+10 1.040 2159.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 8.490000e+10 8.900000e-01 2.503000e+03," &
        " 1.480000e+13 2.600000e-01 4.103000e+03, 2.480000e+15" &
        " -3.300000e-01 6.436000e+03, 2.350000e+25 -3.230000e+00" &
        " 1.316500e+04, 1.020000e+24 -2.670000e+00 1.555200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_95_rev {
  chem_eq = " H + C3H4_P --> H + C3H4_A "
  arrhenius_coeff = 2.4400e+10 1.040 2159.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 8.490000e+10 8.900000e-01 2.503000e+03," &
        " 1.480000e+13 2.600000e-01 4.103000e+03, 2.480000e+15" &
        " -3.300000e-01 6.436000e+03, 2.350000e+25 -3.230000e+00" &
        " 1.316500e+04, 1.020000e+24 -2.670000e+00 1.555200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_96 {
  chem_eq = "H + C3H4_A --> H + C3H4_P"
  arrhenius_coeff = 2.4400e+10 1.040 2159.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.000000e-10 0.000000e+00 0.000000e+00," &
        " 1.000000e-10 0.000000e+00 0.000000e+00, 1.000000e-10" &
        " 0.000000e+00 0.000000e+00, 1.790000e+07 1.980000e+00" &
        " 4.521000e+03, 4.630000e+04 2.620000e+00 4.466000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_96_rev {
  chem_eq = " H + C3H4_P --> H + C3H4_A "
  arrhenius_coeff = 2.4400e+10 1.040 2159.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.000000e-10 0.000000e+00 0.000000e+00," &
        " 1.000000e-10 0.000000e+00 0.000000e+00, 1.000000e-10" &
        " 0.000000e+00 0.000000e+00, 1.790000e+07 1.980000e+00" &
        " 4.521000e+03, 4.630000e+04 2.620000e+00 4.466000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_97 {
  chem_eq = "H + C3H4_A --> C3H5_A"
  arrhenius_coeff = 2.2100e+61 -15.250 20076.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.210000e+61 -1.525000e+01 2.007600e+04," &
        " 1.240000e+52 -1.202000e+01 1.783900e+04, 4.670000e+51" &
        " -1.145000e+01 2.134000e+04, 3.750000e+48 -1.027000e+01" &
        " 2.251100e+04, 4.230000e+43 -8.610000e+00 2.252200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_97_rev {
  chem_eq = " C3H5_A --> H + C3H4_A "
  arrhenius_coeff = 2.2100e+61 -15.250 20076.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.210000e+61 -1.525000e+01 2.007600e+04," &
        " 1.240000e+52 -1.202000e+01 1.783900e+04, 4.670000e+51" &
        " -1.145000e+01 2.134000e+04, 3.750000e+48 -1.027000e+01" &
        " 2.251100e+04, 4.230000e+43 -8.610000e+00 2.252200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_98 {
  chem_eq = "H + C3H4_A --> C3H5_A"
  arrhenius_coeff = 2.2100e+61 -15.250 20076.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.800000e+38 -8.670000e+00 8.035000e+03," &
        " 9.330000e+36 -8.190000e+00 7.462000e+03, 3.320000e+30" &
        " -5.780000e+00 6.913000e+03, 2.290000e+26 -4.320000e+00" &
        " 6.163000e+03, 4.380000e+21 -2.710000e+00 5.187000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_98_rev {
  chem_eq = " C3H5_A --> H + C3H4_A "
  arrhenius_coeff = 2.2100e+61 -15.250 20076.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.800000e+38 -8.670000e+00 8.035000e+03," &
        " 9.330000e+36 -8.190000e+00 7.462000e+03, 3.320000e+30" &
        " -5.780000e+00 6.913000e+03, 2.290000e+26 -4.320000e+00" &
        " 6.163000e+03, 4.380000e+21 -2.710000e+00 5.187000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_99 {
  chem_eq = "H + C3H4_A --> C3H5_S"
  arrhenius_coeff = 1.7600e+10 1.090 3949.00
  press_rxn_param = PLOG "press_coeff: 4.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 3.380000e+49 -1.275000e+01 1.407200e+04, 1.370000e+51" &
        " -1.255000e+01 1.542800e+04, 3.880000e+50 -1.190000e+01" &
        " 1.691500e+04, 2.170000e+49 -1.110000e+01 1.874600e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_99_rev {
  chem_eq = " C3H5_S --> H + C3H4_A "
  arrhenius_coeff = 1.7600e+10 1.090 3949.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 4.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 3.380000e+49 -1.275000e+01 1.407200e+04, 1.370000e+51" &
        " -1.255000e+01 1.542800e+04, 3.880000e+50 -1.190000e+01" &
        " 1.691500e+04, 2.170000e+49 -1.110000e+01 1.874600e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_100 {
  chem_eq = "H + C3H4_A --> C3H5_S"
  arrhenius_coeff = 1.7600e+10 1.090 3949.00
  press_rxn_param = PLOG "press_coeff: 4.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 2.980000e+43 -1.143000e+01 8.046000e+03, 5.750000e+39" &
        " -9.510000e+00 7.458000e+03, 4.330000e+40 -9.600000e+00" &
        " 6.722000e+03, 3.440000e+34 -7.360000e+00 6.150000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_100_rev {
  chem_eq = " C3H5_S --> H + C3H4_A "
  arrhenius_coeff = 1.7600e+10 1.090 3949.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 4.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 2.980000e+43 -1.143000e+01 8.046000e+03, 5.750000e+39" &
        " -9.510000e+00 7.458000e+03, 4.330000e+40 -9.600000e+00" &
        " 6.722000e+03, 3.440000e+34 -7.360000e+00 6.150000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_101 {
  chem_eq = "H + C3H4_A --> C3H5_T"
  arrhenius_coeff = 6.4400e+102 -27.510 51768.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 6.440000e+102 -2.751000e+01 5.176800e+04," &
        " 1.550000e+53 -1.310000e+01 1.447200e+04, 1.900000e+53" &
        " -1.259000e+01 1.672600e+04, 7.950000e+51 -1.182000e+01" &
        " 1.828600e+04, 4.210000e+52 -1.164000e+01 2.226200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_101_rev {
  chem_eq = " C3H5_T --> H + C3H4_A "
  arrhenius_coeff = 6.4400e+102 -27.510 51768.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 6.440000e+102 -2.751000e+01 5.176800e+04," &
        " 1.550000e+53 -1.310000e+01 1.447200e+04, 1.900000e+53" &
        " -1.259000e+01 1.672600e+04, 7.950000e+51 -1.182000e+01" &
        " 1.828600e+04, 4.210000e+52 -1.164000e+01 2.226200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_102 {
  chem_eq = "H + C3H4_A --> C3H5_T"
  arrhenius_coeff = 6.4400e+102 -27.510 51768.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.100000e+54 -1.429000e+01 1.080900e+04," &
        " 9.880000e+44 -1.121000e+01 8.212000e+03, 2.810000e+40" &
        " -9.420000e+00 7.850000e+03, 2.600000e+35 -7.570000e+00" &
        " 7.147000e+03, 9.880000e+29 -5.530000e+00 6.581000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_102_rev {
  chem_eq = " C3H5_T --> H + C3H4_A "
  arrhenius_coeff = 6.4400e+102 -27.510 51768.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.100000e+54 -1.429000e+01 1.080900e+04," &
        " 9.880000e+44 -1.121000e+01 8.212000e+03, 2.810000e+40" &
        " -9.420000e+00 7.850000e+03, 2.600000e+35 -7.570000e+00" &
        " 7.147000e+03, 9.880000e+29 -5.530000e+00 6.581000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_103 {
  chem_eq = "H + C3H4_A --> CH3 + C2H2"
  arrhenius_coeff = 3.7400e+01 3.350 57.80
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.230000e+08 1.530000e+00 4.737000e+03," &
        " 2.720000e+09 1.200000e+00 6.834000e+03, 1.260000e+20" &
        " -1.830000e+00 1.500300e+04, 1.680000e+16 -6.000000e-01" &
        " 1.475400e+04, 1.370000e+17 -7.900000e-01 1.760300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_103_rev {
  chem_eq = " CH3 + C2H2 --> H + C3H4_A "
  arrhenius_coeff = 3.7400e+01 3.350 57.80
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.230000e+08 1.530000e+00 4.737000e+03," &
        " 2.720000e+09 1.200000e+00 6.834000e+03, 1.260000e+20" &
        " -1.830000e+00 1.500300e+04, 1.680000e+16 -6.000000e-01" &
        " 1.475400e+04, 1.370000e+17 -7.900000e-01 1.760300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_104 {
  chem_eq = "H + C3H4_A --> CH3 + C2H2"
  arrhenius_coeff = 3.7400e+01 3.350 57.80
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.000000e-10 0.000000e+00 0.000000e+00," &
        " 1.000000e-10 0.000000e+00 0.000000e+00, 1.230000e+04" &
        " 2.680000e+00 6.335000e+03, 3.310000e+08 1.140000e+00" &
        " 8.886000e+03, 1.280000e+06 1.710000e+00 9.774000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_104_rev {
  chem_eq = " CH3 + C2H2 --> H + C3H4_A "
  arrhenius_coeff = 3.7400e+01 3.350 57.80
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.000000e-10 0.000000e+00 0.000000e+00," &
        " 1.000000e-10 0.000000e+00 0.000000e+00, 1.230000e+04" &
        " 2.680000e+00 6.335000e+03, 3.310000e+08 1.140000e+00" &
        " 8.886000e+03, 1.280000e+06 1.710000e+00 9.774000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_105 {
  chem_eq = "H + C3H4_P --> C3H5_T"
  arrhenius_coeff = 8.8500e+51 -13.040 12325.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 8.850000e+51 -1.304000e+01 1.232500e+04," &
        " 3.170000e+52 -1.269000e+01 1.422600e+04, 2.870000e+53" &
        " -1.251000e+01 1.685300e+04, 9.510000e+51 -1.174000e+01" &
        " 1.833100e+04, 4.510000e+52 -1.158000e+01 2.220700e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_105_rev {
  chem_eq = " C3H5_T --> H + C3H4_P "
  arrhenius_coeff = 8.8500e+51 -13.040 12325.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 8.850000e+51 -1.304000e+01 1.232500e+04," &
        " 3.170000e+52 -1.269000e+01 1.422600e+04, 2.870000e+53" &
        " -1.251000e+01 1.685300e+04, 9.510000e+51 -1.174000e+01" &
        " 1.833100e+04, 4.510000e+52 -1.158000e+01 2.220700e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_106 {
  chem_eq = "H + C3H4_P --> C3H5_T"
  arrhenius_coeff = 8.8500e+51 -13.040 12325.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.970000e+46 -1.191000e+01 7.456000e+03," &
        " 2.590000e+45 -1.123000e+01 8.046000e+03, 6.930000e+39" &
        " -9.110000e+00 7.458000e+03, 6.800000e+34 -7.290000e+00" &
        " 6.722000e+03, 5.650000e+29 -5.390000e+00 6.150000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_106_rev {
  chem_eq = " C3H5_T --> H + C3H4_P "
  arrhenius_coeff = 8.8500e+51 -13.040 12325.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.970000e+46 -1.191000e+01 7.456000e+03," &
        " 2.590000e+45 -1.123000e+01 8.046000e+03, 6.930000e+39" &
        " -9.110000e+00 7.458000e+03, 6.800000e+34 -7.290000e+00" &
        " 6.722000e+03, 5.650000e+29 -5.390000e+00 6.150000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_107 {
  chem_eq = "H + C3H4_P --> C3H5_S"
  arrhenius_coeff = 3.3800e+49 -12.750 14072.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.000000e-10 0.000000e+00 0.000000e+00," &
        " 3.380000e+49 -1.275000e+01 1.407200e+04, 1.370000e+51" &
        " -1.255000e+01 1.542800e+04, 3.880000e+50 -1.190000e+01" &
        " 1.691500e+04, 2.170000e+49 -1.110000e+01 1.874600e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_107_rev {
  chem_eq = " C3H5_S --> H + C3H4_P "
  arrhenius_coeff = 3.3800e+49 -12.750 14072.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.000000e-10 0.000000e+00 0.000000e+00," &
        " 3.380000e+49 -1.275000e+01 1.407200e+04, 1.370000e+51" &
        " -1.255000e+01 1.542800e+04, 3.880000e+50 -1.190000e+01" &
        " 1.691500e+04, 2.170000e+49 -1.110000e+01 1.874600e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_108 {
  chem_eq = "H + C3H4_P --> C3H5_S"
  arrhenius_coeff = 3.3800e+49 -12.750 14072.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.490000e+38 -1.011000e+01 7.458000e+03," &
        " 2.980000e+43 -1.143000e+01 8.736000e+03, 5.750000e+39" &
        " -9.510000e+00 8.772000e+03, 4.330000e+40 -9.600000e+00" &
        " 9.401000e+03, 3.440000e+34 -7.360000e+00 8.558000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_108_rev {
  chem_eq = " C3H5_S --> H + C3H4_P "
  arrhenius_coeff = 3.3800e+49 -12.750 14072.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.490000e+38 -1.011000e+01 7.458000e+03," &
        " 2.980000e+43 -1.143000e+01 8.736000e+03, 5.750000e+39" &
        " -9.510000e+00 8.772000e+03, 4.330000e+40 -9.600000e+00" &
        " 9.401000e+03, 3.440000e+34 -7.360000e+00 8.558000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_109 {
  chem_eq = "H + C3H4_P --> CH3 + C2H2"
  arrhenius_coeff = 2.1200e+10 1.060 3945.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.440000e+10 1.040000e+00 3.980000e+03," &
        " 3.890000e+10 9.890000e-01 4.114000e+03, 3.460000e+12" &
        " 4.420000e-01 5.463000e+03, 1.720000e+14 -1.000000e-02" &
        " 7.134000e+03, 1.900000e+15 -2.900000e-01 8.306000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_109_rev {
  chem_eq = " CH3 + C2H2 --> H + C3H4_P "
  arrhenius_coeff = 2.1200e+10 1.060 3945.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 2.440000e+10 1.040000e+00 3.980000e+03," &
        " 3.890000e+10 9.890000e-01 4.114000e+03, 3.460000e+12" &
        " 4.420000e-01 5.463000e+03, 1.720000e+14 -1.000000e-02" &
        " 7.134000e+03, 1.900000e+15 -2.900000e-01 8.306000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_110 {
  chem_eq = "H + C3H4_P --> C3H5_A"
  arrhenius_coeff = 1.1000e+60 -14.560 28100.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 2.000000e+00 5.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.100000e+60 -1.456000e+01" &
        " 2.810000e+04, 4.910000e+60 -1.437000e+01 3.164400e+04," &
        " 3.040000e+60 -1.419000e+01 3.264200e+04, 9.020000e+59" &
        " -1.389000e+01 3.395300e+04, 2.200000e+59 -1.361000e+01" &
        " 3.490000e+04, 1.600000e+55 -1.207000e+01 3.750000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_110_rev {
  chem_eq = " C3H5_A --> H + C3H4_P "
  arrhenius_coeff = 1.1000e+60 -14.560 28100.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 2.000000e+00 5.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.100000e+60 -1.456000e+01" &
        " 2.810000e+04, 4.910000e+60 -1.437000e+01 3.164400e+04," &
        " 3.040000e+60 -1.419000e+01 3.264200e+04, 9.020000e+59" &
        " -1.389000e+01 3.395300e+04, 2.200000e+59 -1.361000e+01" &
        " 3.490000e+04, 1.600000e+55 -1.207000e+01 3.750000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_111 {
  chem_eq = "C3H5_A --> C3H5_T"
  arrhenius_coeff = 3.9000e+59 -15.420 75400.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 2.000000e+00 5.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.900000e+59 -1.542000e+01" &
        " 7.540000e+04, 7.060000e+56 -1.408000e+01 7.586800e+04," &
        " 4.800000e+55 -1.359000e+01 7.594900e+04, 4.860000e+53" &
        " -1.281000e+01 7.588300e+04, 6.400000e+51 -1.212000e+01" &
        " 7.570000e+04, 2.800000e+43 -9.270000e+00 7.400000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_111_rev {
  chem_eq = " C3H5_T --> C3H5_A "
  arrhenius_coeff = 3.9000e+59 -15.420 75400.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 2.000000e+00 5.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.900000e+59 -1.542000e+01" &
        " 7.540000e+04, 7.060000e+56 -1.408000e+01 7.586800e+04," &
        " 4.800000e+55 -1.359000e+01 7.594900e+04, 4.860000e+53" &
        " -1.281000e+01 7.588300e+04, 6.400000e+51 -1.212000e+01" &
        " 7.570000e+04, 2.800000e+43 -9.270000e+00 7.400000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_112 {
  chem_eq = "C3H5_A --> C3H5_S"
  arrhenius_coeff = 1.3000e+55 -14.530 73800.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 1.300000e+55 -1.453000e+01 7.380000e+04, 5.000000e+51" &
        " -1.302000e+01 7.330000e+04, 9.700000e+48 -1.173000e+01" &
        " 7.370000e+04, 4.860000e+44 -9.840000e+00 7.340000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_112_rev {
  chem_eq = " C3H5_S --> C3H5_A "
  arrhenius_coeff = 1.3000e+55 -14.530 73800.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 1.300000e+55 -1.453000e+01 7.380000e+04, 5.000000e+51" &
        " -1.302000e+01 7.330000e+04, 9.700000e+48 -1.173000e+01" &
        " 7.370000e+04, 4.860000e+44 -9.840000e+00 7.340000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_113 {
  chem_eq = "CH3 + C2H2 --> C3H5_T"
  arrhenius_coeff = 6.8000e+20 -4.160 18000.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 2.000000e+00 5.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 6.800000e+20 -4.160000e+00" &
        " 1.800000e+04, 4.990000e+22 -4.390000e+00 1.885000e+04," &
        " 6.000000e+23 -4.600000e+00 1.957100e+04, 7.310000e+25" &
        " -5.060000e+00 2.115000e+04, 9.300000e+27 -5.550000e+00" &
        " 2.290000e+04, 3.800000e+36 -7.580000e+00 3.130000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_113_rev {
  chem_eq = " C3H5_T --> CH3 + C2H2 "
  arrhenius_coeff = 6.8000e+20 -4.160 18000.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 2.000000e+00 5.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 6.800000e+20 -4.160000e+00" &
        " 1.800000e+04, 4.990000e+22 -4.390000e+00 1.885000e+04," &
        " 6.000000e+23 -4.600000e+00 1.957100e+04, 7.310000e+25" &
        " -5.060000e+00 2.115000e+04, 9.300000e+27 -5.550000e+00" &
        " 2.290000e+04, 3.800000e+36 -7.580000e+00 3.130000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_114 {
  chem_eq = "C3H5_T --> C3H5_S"
  arrhenius_coeff = 1.6000e+44 -12.160 52200.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 1.600000e+44 -1.216000e+01 5.220000e+04, 1.500000e+48" &
        " -1.271000e+01 5.390000e+04, 5.100000e+52 -1.337000e+01" &
        " 5.720000e+04, 5.800000e+51 -1.243000e+01 5.920000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_114_rev {
  chem_eq = " C3H5_S --> C3H5_T "
  arrhenius_coeff = 1.6000e+44 -12.160 52200.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 1.600000e+44 -1.216000e+01 5.220000e+04, 1.500000e+48" &
        " -1.271000e+01 5.390000e+04, 5.100000e+52 -1.337000e+01" &
        " 5.720000e+04, 5.800000e+51 -1.243000e+01 5.920000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_115 {
  chem_eq = "CH3 + C2H2 --> C3H5_A"
  arrhenius_coeff = 8.2000e+53 -13.320 33200.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 2.000000e+00 5.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 8.200000e+53 -1.332000e+01" &
        " 3.320000e+04, 2.680000e+53 -1.282000e+01 3.573000e+04," &
        " 3.640000e+52 -1.246000e+01 3.612700e+04, 1.040000e+51" &
        " -1.189000e+01 3.647600e+04, 4.400000e+49 -1.140000e+01" &
        " 3.670000e+04, 3.800000e+44 -9.630000e+00 3.760000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_115_rev {
  chem_eq = " C3H5_A --> CH3 + C2H2 "
  arrhenius_coeff = 8.2000e+53 -13.320 33200.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 2.000000e+00 5.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 8.200000e+53 -1.332000e+01" &
        " 3.320000e+04, 2.680000e+53 -1.282000e+01 3.573000e+04," &
        " 3.640000e+52 -1.246000e+01 3.612700e+04, 1.040000e+51" &
        " -1.189000e+01 3.647600e+04, 4.400000e+49 -1.140000e+01" &
        " 3.670000e+04, 3.800000e+44 -9.630000e+00 3.760000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_116 {
  chem_eq = "CH3 + C2H2 --> C3H5_S"
  arrhenius_coeff = 1.7800e+42 -10.400 13647.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.780000e+42 -1.040000e+01 1.364700e+04," &
        " 1.520000e+44 -1.073000e+01 1.525600e+04, 1.190000e+44" &
        " -1.019000e+01 1.872800e+04, 6.020000e+43 -9.740000e+00" &
        " 2.056100e+04, 1.420000e+42 -8.910000e+00 2.223500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_116_rev {
  chem_eq = " C3H5_S --> CH3 + C2H2 "
  arrhenius_coeff = 1.7800e+42 -10.400 13647.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.780000e+42 -1.040000e+01 1.364700e+04," &
        " 1.520000e+44 -1.073000e+01 1.525600e+04, 1.190000e+44" &
        " -1.019000e+01 1.872800e+04, 6.020000e+43 -9.740000e+00" &
        " 2.056100e+04, 1.420000e+42 -8.910000e+00 2.223500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_117 {
  chem_eq = "CH3 + C2H2 --> C3H5_S"
  arrhenius_coeff = 1.7800e+42 -10.400 13647.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.000000e-10 0.000000e+00 0.000000e+00," &
        " 1.000000e-10 0.000000e+00 0.000000e+00, 8.490000e+35" &
        " -8.430000e+00 1.235600e+04, 3.040000e+32 -7.010000e+00" &
        " 1.235700e+04, 1.690000e+27 -5.070000e+00 1.169000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_117_rev {
  chem_eq = " C3H5_S --> CH3 + C2H2 "
  arrhenius_coeff = 1.7800e+42 -10.400 13647.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 3.900000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.000000e-10 0.000000e+00 0.000000e+00," &
        " 1.000000e-10 0.000000e+00 0.000000e+00, 8.490000e+35" &
        " -8.430000e+00 1.235600e+04, 3.040000e+32 -7.010000e+00" &
        " 1.235700e+04, 1.690000e+27 -5.070000e+00 1.169000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_118 {
  chem_eq = "C2H + C3H4_A --> C2H2 + C3H3"
  arrhenius_coeff = 1.0000e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_118_rev {
  chem_eq = " C2H2 + C3H3 --> C2H + C3H4_A "
  arrhenius_coeff = 1.0000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_119 {
  chem_eq = "C2H5 + C2H --> CH3 + C3H3"
  arrhenius_coeff = 1.8100e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_119_rev {
  chem_eq = " CH3 + C3H3 --> C2H5 + C2H "
  arrhenius_coeff = 1.8100e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_120 {
  chem_eq = "NC4H10 --> 2.0*C2H5"
  arrhenius_coeff = 1.3550e+37 -6.036 92929.00
  press_rxn_param = Troe_falloff "arrhenius_press: 4.72e+18" &
        " 0.000 49578.0 press_coeff: 0.07998 1.000e-20 3.243e+04" &
        " 4858."
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_120_rev {
  chem_eq = " 2.0*C2H5 --> NC4H10 "
  arrhenius_coeff = 1.3550e+37 -6.036 92929.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 4.72e+18" &
        " 0.000 49578.0 press_coeff: 0.07998 1.000e-20 3.243e+04" &
        " 4858."
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_121 {
  chem_eq = "NC4H10 --> CH3 + NC3H7"
  arrhenius_coeff = 6.6000e+52 -10.626 100330.00
  press_rxn_param = Troe_falloff "arrhenius_press: 5.34e+17" &
        " 0.000 42959.0 press_coeff: 0.09502 1.000e-20 5348. 4326."
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_121_rev {
  chem_eq = " CH3 + NC3H7 --> NC4H10 "
  arrhenius_coeff = 6.6000e+52 -10.626 100330.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 5.34e+17" &
        " 0.000 42959.0 press_coeff: 0.09502 1.000e-20 5348. 4326."
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_122 {
  chem_eq = "NC4H10 --> H + PC4H9"
  arrhenius_coeff = 4.8000e+40 -7.060 115302.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 4.450000e+90 -2.191000e+01 1.405640e+05," &
        " 4.630000e+76 -1.764000e+01 1.346690e+05, 4.940000e+58" &
        " -1.232000e+01 1.254350e+05, 4.800000e+40 -7.060000e+00" &
        " 1.153020e+05, 1.490000e+27 -3.150000e+00 1.073230e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_122_rev {
  chem_eq = " H + PC4H9 --> NC4H10 "
  arrhenius_coeff = 4.8000e+40 -7.060 115302.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 4.450000e+90 -2.191000e+01 1.405640e+05," &
        " 4.630000e+76 -1.764000e+01 1.346690e+05, 4.940000e+58" &
        " -1.232000e+01 1.254350e+05, 4.800000e+40 -7.060000e+00" &
        " 1.153020e+05, 1.490000e+27 -3.150000e+00 1.073230e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_123 {
  chem_eq = "NC4H10 --> H + SC4H9"
  arrhenius_coeff = 8.5200e+38 -6.580 110556.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 3.100000e+88 -2.124000e+01 1.363550e+05," &
        " 4.340000e+73 -1.676000e+01 1.295900e+05, 7.390000e+55" &
        " -1.152000e+01 1.201990e+05, 8.520000e+38 -6.580000e+00" &
        " 1.105560e+05, 5.400000e+26 -3.050000e+00 1.033130e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_123_rev {
  chem_eq = " H + SC4H9 --> NC4H10 "
  arrhenius_coeff = 8.5200e+38 -6.580 110556.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 3.100000e+88 -2.124000e+01 1.363550e+05," &
        " 4.340000e+73 -1.676000e+01 1.295900e+05, 7.390000e+55" &
        " -1.152000e+01 1.201990e+05, 8.520000e+38 -6.580000e+00" &
        " 1.105560e+05, 5.400000e+26 -3.050000e+00 1.033130e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_124 {
  chem_eq = "H + NC4H10 --> H2 + PC4H9"
  arrhenius_coeff = 1.7500e+05 2.690 6450.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_124_rev {
  chem_eq = " H2 + PC4H9 --> H + NC4H10 "
  arrhenius_coeff = 1.7500e+05 2.690 6450.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_125 {
  chem_eq = "H + NC4H10 --> H2 + SC4H9"
  arrhenius_coeff = 1.3000e+06 2.400 4471.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_125_rev {
  chem_eq = " H2 + SC4H9 --> H + NC4H10 "
  arrhenius_coeff = 1.3000e+06 2.400 4471.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_126 {
  chem_eq = "IC4H10 --> CH3 + IC3H7"
  arrhenius_coeff = 2.5200e+31 -4.102 91495.00
  press_rxn_param = Troe_falloff "arrhenius_press: 2.41e+19" &
        " 0.000 52576.0 press_coeff: 0.3662 815.3 60.79 1.000e+20"
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_126_rev {
  chem_eq = " CH3 + IC3H7 --> IC4H10 "
  arrhenius_coeff = 2.5200e+31 -4.102 91495.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 2.41e+19" &
        " 0.000 52576.0 press_coeff: 0.3662 815.3 60.79 1.000e+20"
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_127 {
  chem_eq = "IC4H10 --> H + IC4H9"
  arrhenius_coeff = 9.8500e+95 -23.110 147600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_127_rev {
  chem_eq = " H + IC4H9 --> IC4H10 "
  arrhenius_coeff = 9.8500e+95 -23.110 147600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_128 {
  chem_eq = "CH3 + IC4H10 --> CH4 + IC4H9"
  arrhenius_coeff = 1.3600e+00 3.650 7154.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_128_rev {
  chem_eq = " CH4 + IC4H9 --> CH3 + IC4H10 "
  arrhenius_coeff = 1.3600e+00 3.650 7154.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_129 {
  chem_eq = "H + IC4H10 --> H2 + IC4H9"
  arrhenius_coeff = 1.8100e+06 2.540 6756.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_129_rev {
  chem_eq = " H2 + IC4H9 --> H + IC4H10 "
  arrhenius_coeff = 1.8100e+06 2.540 6756.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_130 {
  chem_eq = "C2H5 + IC4H10 --> C2H6 + IC4H9"
  arrhenius_coeff = 1.5100e+12 0.000 10400.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_130_rev {
  chem_eq = " C2H6 + IC4H9 --> C2H5 + IC4H10 "
  arrhenius_coeff = 1.5100e+12 0.000 10400.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_131 {
  chem_eq = "IC4H9 --> CH3 + C3H6"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 3.150000e+41 -9.500000e+00 3.348600e+04, 6.750000e+44" &
        " -1.007000e+01 3.720900e+04, 7.790000e+44 -9.700000e+00" &
        " 3.975100e+04, 3.610000e+39 -7.780000e+00 3.958300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_131_rev {
  chem_eq = " CH3 + C3H6 --> IC4H9 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 3.150000e+41 -9.500000e+00 3.348600e+04, 6.750000e+44" &
        " -1.007000e+01 3.720900e+04, 7.790000e+44 -9.700000e+00" &
        " 3.975100e+04, 3.610000e+39 -7.780000e+00 3.958300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_132 {
  chem_eq = "IC4H8 --> CH3 + C3H5_T"
  arrhenius_coeff = 1.4200e+93 -22.790 133825.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 3.500000e+00 1.000000e+01 3.500000e+01" &
        " 1.000000e+02 arrhenius_press: 1.260000e+94 -2.299000e+01" &
        " 1.340240e+05, 6.760000e+93 -2.251000e+01 1.379330e+05," &
        " 3.140000e+90 -2.137000e+01 1.378660e+05, 9.200000e+85" &
        " -1.994000e+01 1.364980e+05, 6.050000e+78 -1.776000e+01" &
        " 1.334370e+05, 4.870000e+71 -1.565000e+01 1.299190e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_132_rev {
  chem_eq = " CH3 + C3H5_T --> IC4H8 "
  arrhenius_coeff = 1.4200e+93 -22.790 133825.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 3.500000e+00 1.000000e+01 3.500000e+01" &
        " 1.000000e+02 arrhenius_press: 1.260000e+94 -2.299000e+01" &
        " 1.340240e+05, 6.760000e+93 -2.251000e+01 1.379330e+05," &
        " 3.140000e+90 -2.137000e+01 1.378660e+05, 9.200000e+85" &
        " -1.994000e+01 1.364980e+05, 6.050000e+78 -1.776000e+01" &
        " 1.334370e+05, 4.870000e+71 -1.565000e+01 1.299190e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_133 {
  chem_eq = "IC4H8 --> H + IC4H7"
  arrhenius_coeff = 4.6600e+92 -22.450 129059.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 3.500000e+00 1.000000e+01 3.500000e+01" &
        " 1.000000e+02 arrhenius_press: 7.510000e+95 -2.338000e+01" &
        " 1.292140e+05, 3.590000e+88 -2.099000e+01 1.278130e+05," &
        " 2.960000e+82 -1.912000e+01 1.254560e+05, 2.130000e+76" &
        " -1.727000e+01 1.226290e+05, 1.130000e+68 -1.482000e+01" &
        " 1.184160e+05, 4.730000e+60 -1.266000e+01 1.144040e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_133_rev {
  chem_eq = " H + IC4H7 --> IC4H8 "
  arrhenius_coeff = 4.6600e+92 -22.450 129059.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 3.500000e+00 1.000000e+01 3.500000e+01" &
        " 1.000000e+02 arrhenius_press: 7.510000e+95 -2.338000e+01" &
        " 1.292140e+05, 3.590000e+88 -2.099000e+01 1.278130e+05," &
        " 2.960000e+82 -1.912000e+01 1.254560e+05, 2.130000e+76" &
        " -1.727000e+01 1.226290e+05, 1.130000e+68 -1.482000e+01" &
        " 1.184160e+05, 4.730000e+60 -1.266000e+01 1.144040e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_134 {
  chem_eq = "C3H5_A + IC4H8 --> C3H6 + IC4H7"
  arrhenius_coeff = 7.9400e+11 0.000 20500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_134_rev {
  chem_eq = " C3H6 + IC4H7 --> C3H5_A + IC4H8 "
  arrhenius_coeff = 7.9400e+11 0.000 20500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_135 {
  chem_eq = "C3H5_S + IC4H8 --> C3H6 + IC4H7"
  arrhenius_coeff = 7.9400e+11 0.000 20500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_135_rev {
  chem_eq = " C3H6 + IC4H7 --> C3H5_S + IC4H8 "
  arrhenius_coeff = 7.9400e+11 0.000 20500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_136 {
  chem_eq = "C3H5_T + IC4H8 --> C3H6 + IC4H7"
  arrhenius_coeff = 7.9400e+11 0.000 20500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_136_rev {
  chem_eq = " C3H6 + IC4H7 --> C3H5_T + IC4H8 "
  arrhenius_coeff = 7.9400e+11 0.000 20500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_137 {
  chem_eq = "CH3 + C3H4_A --> IC4H7"
  arrhenius_coeff = 4.0210e+04 2.500 8847.50
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_137_rev {
  chem_eq = " IC4H7 --> CH3 + C3H4_A "
  arrhenius_coeff = 4.0210e+04 2.500 8847.50
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_138 {
  chem_eq = "H + IC4H8 --> IC4H9"
  arrhenius_coeff = 1.0000e+00 1.000 0.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 arrhenius_press:" &
        " 7.990000e+81 -2.316100e+01 2.223900e+04, 4.240000e+68" &
        " -1.842700e+01 1.966500e+04, 1.040000e+49 -1.150000e+01" &
        " 1.535900e+04, 6.200000e+41 -8.892000e+00 1.463700e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_138_rev {
  chem_eq = " IC4H9 --> H + IC4H8 "
  arrhenius_coeff = 1.0000e+00 1.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 arrhenius_press:" &
        " 7.990000e+81 -2.316100e+01 2.223900e+04, 4.240000e+68" &
        " -1.842700e+01 1.966500e+04, 1.040000e+49 -1.150000e+01" &
        " 1.535900e+04, 6.200000e+41 -8.892000e+00 1.463700e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_139 {
  chem_eq = "H + IC4H8 --> IC4H9"
  arrhenius_coeff = 1.0000e+00 1.000 0.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.850000e+26 -5.830000e+00 3.865800e+03," &
        " 2.820000e+30 -6.490000e+00 5.470800e+03, 3.780000e+28" &
        " -5.570000e+00 5.625100e+03, 1.460000e+25 -4.280000e+00" &
        " 5.247800e+03, 4.220000e+27 -4.390000e+00 9.345800e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_139_rev {
  chem_eq = " IC4H9 --> H + IC4H8 "
  arrhenius_coeff = 1.0000e+00 1.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.850000e+26 -5.830000e+00 3.865800e+03," &
        " 2.820000e+30 -6.490000e+00 5.470800e+03, 3.780000e+28" &
        " -5.570000e+00 5.625100e+03, 1.460000e+25 -4.280000e+00" &
        " 5.247800e+03, 4.220000e+27 -4.390000e+00 9.345800e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_140 {
  chem_eq = "H + IC4H8 --> CH3 + C3H6"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 5.130000e+08 1.350000e+00 2.542000e+03," &
        " 2.630000e+10 8.700000e-01 3.599600e+03, 8.900000e+11" &
        " 4.700000e-01 5.431100e+03, 3.080000e+22 -2.600000e+00" &
        " 1.289800e+04, 4.400000e+22 -2.420000e+00 1.650000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_140_rev {
  chem_eq = " CH3 + C3H6 --> H + IC4H8 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 5.130000e+08 1.350000e+00 2.542000e+03," &
        " 2.630000e+10 8.700000e-01 3.599600e+03, 8.900000e+11" &
        " 4.700000e-01 5.431100e+03, 3.080000e+22 -2.600000e+00" &
        " 1.289800e+04, 4.400000e+22 -2.420000e+00 1.650000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_141 {
  chem_eq = "H + IC4H8 --> CH3 + C3H6"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 7.700000e+02 1.350000e+00 2.542000e+03," &
        " 3.940000e+04 8.700000e-01 3.599600e+03, 1.340000e+06" &
        " 4.700000e-01 5.431100e+03, 4.130000e+04 2.520000e+00" &
        " 3.679100e+03, 8.370000e+02 2.910000e+00 3.980900e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_141_rev {
  chem_eq = " CH3 + C3H6 --> H + IC4H8 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.300000e-03" &
        " 4.000000e-02 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 7.700000e+02 1.350000e+00 2.542000e+03," &
        " 3.940000e+04 8.700000e-01 3.599600e+03, 1.340000e+06" &
        " 4.700000e-01 5.431100e+03, 4.130000e+04 2.520000e+00" &
        " 3.679100e+03, 8.370000e+02 2.910000e+00 3.980900e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_142 {
  chem_eq = "CH3 + C3H5_A --> C4H8_1"
  arrhenius_coeff = 6.0000e+14 -0.320 -262.30
  press_rxn_param = Troe_falloff "arrhenius_press: 3.91e+60" &
        " -12.810 6250.0 press_coeff: 0.1040 1606. 6.000e+04 6118."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_142_rev {
  chem_eq = " C4H8_1 --> CH3 + C3H5_A "
  arrhenius_coeff = 6.0000e+14 -0.320 -262.30
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 3.91e+60" &
        " -12.810 6250.0 press_coeff: 0.1040 1606. 6.000e+04 6118."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_143 {
  chem_eq = "C2H5 + C2H3 --> C4H8_1"
  arrhenius_coeff = 1.5000e+13 0.000 0.00
  press_rxn_param = Troe_falloff "arrhenius_press: 1.55e+56" &
        " -11.790 8984.5 press_coeff: 0.1980 2278. 6.000e+04 5723."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_143_rev {
  chem_eq = " C4H8_1 --> C2H5 + C2H3 "
  arrhenius_coeff = 1.5000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 1.55e+56" &
        " -11.790 8984.5 press_coeff: 0.1980 2278. 6.000e+04 5723."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_144 {
  chem_eq = "H + C4H71_4 --> C4H8_1"
  arrhenius_coeff = 3.6000e+13 0.000 0.00
  press_rxn_param = Troe_falloff "arrhenius_press: 3.01e+48" &
        " -9.320 5833.6 press_coeff: 0.4980 1314. 1314. 5.000e+04"
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_144_rev {
  chem_eq = " C4H8_1 --> H + C4H71_4 "
  arrhenius_coeff = 3.6000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 3.01e+48" &
        " -9.320 5833.6 press_coeff: 0.4980 1314. 1314. 5.000e+04"
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_145 {
  chem_eq = "H + C4H71_3 --> C4H8_1"
  arrhenius_coeff = 5.0000e+13 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_145_rev {
  chem_eq = " C4H8_1 --> H + C4H71_3 "
  arrhenius_coeff = 5.0000e+13 0.000 5000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_146 {
  chem_eq = "H + C4H71_3 --> C4H8_2"
  arrhenius_coeff = 5.0000e+13 0.000 0.00
  press_rxn_param = Troe_falloff "arrhenius_press: 1.33e+60" &
        " -12.000 5967.8 press_coeff: 0.02000 1097. 1.097e+04 6860."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_146_rev {
  chem_eq = " C4H8_2 --> H + C4H71_3 "
  arrhenius_coeff = 5.0000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 1.33e+60" &
        " -12.000 5967.8 press_coeff: 0.02000 1097. 1.097e+04 6860."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_147 {
  chem_eq = "CH3 + C3H5_S --> C4H8_2"
  arrhenius_coeff = 5.0000e+13 0.000 0.00
  press_rxn_param = Troe_falloff "arrhenius_press: 8.54e+58" &
        " -11.940 9769.8 press_coeff: 0.1750 1341. 6.000e+04" &
        " 1.014e+04"
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00" &
        " C2H2:3.00 C2H4:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_147_rev {
  chem_eq = " C4H8_2 --> CH3 + C3H5_S "
  arrhenius_coeff = 5.0000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 8.54e+58" &
        " -11.940 9769.8 press_coeff: 0.1750 1341. 6.000e+04" &
        " 1.014e+04"
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00" &
        " C2H2:3.00 C2H4:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_148 {
  chem_eq = "C4H8_2 --> CH3 + C3H5_A"
  arrhenius_coeff = 7.5000e+65 -15.600 97300.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_148_rev {
  chem_eq = " CH3 + C3H5_A --> C4H8_2 "
  arrhenius_coeff = 7.5000e+65 -15.600 97300.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_149 {
  chem_eq = "C4H8_2 --> H + C4H71_3"
  arrhenius_coeff = 4.6000e+84 -20.030 132787.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_149_rev {
  chem_eq = " H + C4H71_3 --> C4H8_2 "
  arrhenius_coeff = 4.6000e+84 -20.030 132787.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_150 {
  chem_eq = "H + C4H8_2 --> H2 + C4H71_3"
  arrhenius_coeff = 5.6200e+02 3.500 1627.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_150_rev {
  chem_eq = " H2 + C4H71_3 --> H + C4H8_2 "
  arrhenius_coeff = 5.6200e+02 3.500 1627.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_151 {
  chem_eq = "CH3 + C4H8_2 --> CH4 + C4H71_3"
  arrhenius_coeff = 7.1400e+00 3.570 7642.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_151_rev {
  chem_eq = " CH4 + C4H71_3 --> CH3 + C4H8_2 "
  arrhenius_coeff = 7.1400e+00 3.570 7642.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_152 {
  chem_eq = "C2H5 + C4H71_3 --> C2H4 + C4H8_1"
  arrhenius_coeff = 2.5900e+12 0.000 -131.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_152_rev {
  chem_eq = " C2H4 + C4H8_1 --> C2H5 + C4H71_3 "
  arrhenius_coeff = 2.5900e+12 0.000 -131.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_153 {
  chem_eq = "C4H71_3 --> H + C4H6"
  arrhenius_coeff = 8.5300e+07 1.950 47490.11
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_153_rev {
  chem_eq = " H + C4H6 --> C4H71_3 "
  arrhenius_coeff = 8.5300e+07 1.950 47490.11
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_154 {
  chem_eq = "C4H71_4 --> C2H4 + C2H3"
  arrhenius_coeff = 2.8400e+10 0.990 38998.80
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_154_rev {
  chem_eq = " C2H4 + C2H3 --> C4H71_4 "
  arrhenius_coeff = 2.8400e+10 0.990 38998.80
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_155 {
  chem_eq = "C4H71_4 --> H + C4H6"
  arrhenius_coeff = 1.3200e+05 2.280 33245.86
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_155_rev {
  chem_eq = " H + C4H6 --> C4H71_4 "
  arrhenius_coeff = 1.3200e+05 2.280 33245.86
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_156 {
  chem_eq = "C4H71_3 --> C4H71_4"
  arrhenius_coeff = 5.6200e-12 7.190 36200.82
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_156_rev {
  chem_eq = " C4H71_4 --> C4H71_3 "
  arrhenius_coeff = 5.6200e-12 7.190 36200.82
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_157 {
  chem_eq = "H + C4H8_1 --> C2H5 + C2H4"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 2.550000e+06 1.930000e+00" &
        " 5.564000e+03, 5.560000e+06 1.830000e+00 5.802000e+03," &
        " 1.210000e+09 1.180000e+00 7.472000e+03, 9.470000e+16" &
        " -1.030000e+00 1.341300e+04, 4.500000e+28 -4.240000e+00" &
        " 2.361800e+04, 7.020000e+32 -5.220000e+00 3.175400e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_157_rev {
  chem_eq = " C2H5 + C2H4 --> H + C4H8_1 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 2.550000e+06 1.930000e+00" &
        " 5.564000e+03, 5.560000e+06 1.830000e+00 5.802000e+03," &
        " 1.210000e+09 1.180000e+00 7.472000e+03, 9.470000e+16" &
        " -1.030000e+00 1.341300e+04, 4.500000e+28 -4.240000e+00" &
        " 2.361800e+04, 7.020000e+32 -5.220000e+00 3.175400e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_158 {
  chem_eq = "H + C4H8_1 --> C2H5 + C2H4"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.450000e+07 1.810000e+00" &
        " 2.263000e+03, 8.060000e+07 1.710000e+00 2.522000e+03," &
        " 1.180000e+10 1.100000e+00 4.077000e+03, 6.020000e+15" &
        " -4.900000e-01 8.452000e+03, 7.580000e+21 -2.140000e+00" &
        " 1.424500e+04, 2.290000e+21 -1.870000e+00 1.724300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_158_rev {
  chem_eq = " C2H5 + C2H4 --> H + C4H8_1 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.450000e+07 1.810000e+00" &
        " 2.263000e+03, 8.060000e+07 1.710000e+00 2.522000e+03," &
        " 1.180000e+10 1.100000e+00 4.077000e+03, 6.020000e+15" &
        " -4.900000e-01 8.452000e+03, 7.580000e+21 -2.140000e+00" &
        " 1.424500e+04, 2.290000e+21 -1.870000e+00 1.724300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_159 {
  chem_eq = "H + C4H8_1 --> CH3 + C3H6"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 7.830000e+09 1.170000e+00" &
        " 1.442000e+03, 3.390000e+10 1.000000e+00 1.895000e+03," &
        " 3.700000e+13 1.400000e-01 4.127000e+03, 4.570000e+19" &
        " -1.540000e+00 9.061000e+03, 8.570000e+23 -2.660000e+00" &
        " 1.414000e+04, 1.320000e+20 -1.460000e+00 1.538300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_159_rev {
  chem_eq = " CH3 + C3H6 --> H + C4H8_1 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 7.830000e+09 1.170000e+00" &
        " 1.442000e+03, 3.390000e+10 1.000000e+00 1.895000e+03," &
        " 3.700000e+13 1.400000e-01 4.127000e+03, 4.570000e+19" &
        " -1.540000e+00 9.061000e+03, 8.570000e+23 -2.660000e+00" &
        " 1.414000e+04, 1.320000e+20 -1.460000e+00 1.538300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_160 {
  chem_eq = "H + C4H8_1 --> CH3 + C3H6"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.800000e+06 1.760000e+00" &
        " 5.900000e+03, 3.460000e+06 1.680000e+00 6.100000e+03," &
        " 4.020000e+08 1.100000e+00 7.574000e+03, 1.210000e+16" &
        " -9.900000e-01 1.317500e+04, 7.140000e+27 -4.230000e+00" &
        " 2.331900e+04, 1.000000e+33 -5.490000e+00 3.192200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_160_rev {
  chem_eq = " CH3 + C3H6 --> H + C4H8_1 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.800000e+06 1.760000e+00" &
        " 5.900000e+03, 3.460000e+06 1.680000e+00 6.100000e+03," &
        " 4.020000e+08 1.100000e+00 7.574000e+03, 1.210000e+16" &
        " -9.900000e-01 1.317500e+04, 7.140000e+27 -4.230000e+00" &
        " 2.331900e+04, 1.000000e+33 -5.490000e+00 3.192200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_161 {
  chem_eq = "H + C4H8_1 --> PC4H9"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.350000e+15 -2.810000e+00" &
        " 1.570000e+03, 5.200000e+16 -2.970000e+00 1.992000e+03," &
        " 1.910000e+21 -3.970000e+00 4.636000e+03, 1.900000e+31" &
        " -6.460000e+00 1.196800e+04, 2.100000e+40 -8.600000e+00" &
        " 2.105800e+04, 1.440000e+37 -7.210000e+00 2.489600e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_161_rev {
  chem_eq = " PC4H9 --> H + C4H8_1 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.350000e+15 -2.810000e+00" &
        " 1.570000e+03, 5.200000e+16 -2.970000e+00 1.992000e+03," &
        " 1.910000e+21 -3.970000e+00 4.636000e+03, 1.900000e+31" &
        " -6.460000e+00 1.196800e+04, 2.100000e+40 -8.600000e+00" &
        " 2.105800e+04, 1.440000e+37 -7.210000e+00 2.489600e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_162 {
  chem_eq = "H + C4H8_1 --> PC4H9"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 4.330000e+20 -4.160000e+00" &
        " -2.630000e+02, 1.780000e+22 -4.330000e+00 1.860000e+02," &
        " 1.980000e+26 -5.180000e+00 2.518000e+03, 3.780000e+32" &
        " -6.630000e+00 7.265000e+03, 8.790000e+34 -6.910000e+00" &
        " 1.095200e+04, 7.800000e+28 -4.790000e+00 1.035500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_162_rev {
  chem_eq = " PC4H9 --> H + C4H8_1 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 4.330000e+20 -4.160000e+00" &
        " -2.630000e+02, 1.780000e+22 -4.330000e+00 1.860000e+02," &
        " 1.980000e+26 -5.180000e+00 2.518000e+03, 3.780000e+32" &
        " -6.630000e+00 7.265000e+03, 8.790000e+34 -6.910000e+00" &
        " 1.095200e+04, 7.800000e+28 -4.790000e+00 1.035500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_163 {
  chem_eq = "H + C4H8_1 --> SC4H9"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 4.070000e+22 -4.510000e+00" &
        " -7.710000e+02, 3.900000e+24 -4.780000e+00 -3.400000e+01," &
        " 2.030000e+29 -5.810000e+00 2.970000e+03, 3.530000e+34" &
        " -6.950000e+00 7.525000e+03, 1.190000e+34 -6.420000e+00" &
        " 9.810000e+03, 1.370000e+26 -3.790000e+00 8.012000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_163_rev {
  chem_eq = " SC4H9 --> H + C4H8_1 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 4.070000e+22 -4.510000e+00" &
        " -7.710000e+02, 3.900000e+24 -4.780000e+00 -3.400000e+01," &
        " 2.030000e+29 -5.810000e+00 2.970000e+03, 3.530000e+34" &
        " -6.950000e+00 7.525000e+03, 1.190000e+34 -6.420000e+00" &
        " 9.810000e+03, 1.370000e+26 -3.790000e+00 8.012000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_164 {
  chem_eq = "H + C4H8_1 --> SC4H9"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.520000e+12 -2.150000e+00" &
        " 1.466000e+03, 1.020000e+14 -2.280000e+00 1.799000e+03," &
        " 1.160000e+18 -3.130000e+00 4.049000e+03, 5.220000e+27" &
        " -5.530000e+00 1.096300e+04, 4.330000e+37 -7.920000e+00" &
        " 2.035400e+04, 2.220000e+36 -7.060000e+00 2.520300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_164_rev {
  chem_eq = " SC4H9 --> H + C4H8_1 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.520000e+12 -2.150000e+00" &
        " 1.466000e+03, 1.020000e+14 -2.280000e+00 1.799000e+03," &
        " 1.160000e+18 -3.130000e+00 4.049000e+03, 5.220000e+27" &
        " -5.530000e+00 1.096300e+04, 4.330000e+37 -7.920000e+00" &
        " 2.035400e+04, 2.220000e+36 -7.060000e+00 2.520300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_165 {
  chem_eq = "H + C4H8_2 --> C2H5 + C2H4"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 8.960000e+06 1.860000e+00" &
        " 6.209000e+03, 1.920000e+07 1.770000e+00 6.443000e+03," &
        " 3.970000e+09 1.110000e+00 8.097000e+03, 3.010000e+17" &
        " -1.090000e+00 1.402300e+04, 1.880000e+29 -4.330000e+00" &
        " 2.429700e+04, 5.150000e+33 -5.390000e+00 3.260100e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_165_rev {
  chem_eq = " C2H5 + C2H4 --> H + C4H8_2 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 8.960000e+06 1.860000e+00" &
        " 6.209000e+03, 1.920000e+07 1.770000e+00 6.443000e+03," &
        " 3.970000e+09 1.110000e+00 8.097000e+03, 3.010000e+17" &
        " -1.090000e+00 1.402300e+04, 1.880000e+29 -4.330000e+00" &
        " 2.429700e+04, 5.150000e+33 -5.390000e+00 3.260100e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_166 {
  chem_eq = "H + C4H8_2 --> CH3 + C3H6"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 6.390000e+09 1.290000e+00" &
        " 1.834000e+03, 2.600000e+10 1.120000e+00 2.267000e+03," &
        " 2.480000e+13 2.900000e-01 4.456000e+03, 2.910000e+19" &
        " -1.390000e+00 9.365000e+03, 6.130000e+23 -2.530000e+00" &
        " 1.446300e+04, 1.230000e+20 -1.350000e+00 1.576200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_166_rev {
  chem_eq = " CH3 + C3H6 --> H + C4H8_2 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 6.390000e+09 1.290000e+00" &
        " 1.834000e+03, 2.600000e+10 1.120000e+00 2.267000e+03," &
        " 2.480000e+13 2.900000e-01 4.456000e+03, 2.910000e+19" &
        " -1.390000e+00 9.365000e+03, 6.130000e+23 -2.530000e+00" &
        " 1.446300e+04, 1.230000e+20 -1.350000e+00 1.576200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_167 {
  chem_eq = "H + C4H8_2 --> PC4H9"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.900000e+14 -2.550000e+00" &
        " 1.729000e+03, 1.410000e+16 -2.710000e+00 2.133000e+03," &
        " 4.310000e+20 -3.690000e+00 4.719000e+03, 4.030000e+30" &
        " -6.170000e+00 1.202000e+04, 5.190000e+39 -8.330000e+00" &
        " 2.113700e+04, 5.170000e+36 -6.980000e+00 2.506300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_167_rev {
  chem_eq = " PC4H9 --> H + C4H8_2 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.900000e+14 -2.550000e+00" &
        " 1.729000e+03, 1.410000e+16 -2.710000e+00 2.133000e+03," &
        " 4.310000e+20 -3.690000e+00 4.719000e+03, 4.030000e+30" &
        " -6.170000e+00 1.202000e+04, 5.190000e+39 -8.330000e+00" &
        " 2.113700e+04, 5.170000e+36 -6.980000e+00 2.506300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_168 {
  chem_eq = "H + C4H8_2 --> SC4H9"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 8.340000e+21 -4.210000e+00" &
        " -6.020000e+02, 6.790000e+23 -4.460000e+00 8.200000e+01," &
        " 2.850000e+28 -5.470000e+00 3.003000e+03, 5.450000e+33" &
        " -6.610000e+00 7.559000e+03, 2.330000e+33 -6.110000e+00" &
        " 9.893000e+03, 3.270000e+25 -3.510000e+00 8.145000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_168_rev {
  chem_eq = " SC4H9 --> H + C4H8_2 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 8.340000e+21 -4.210000e+00" &
        " -6.020000e+02, 6.790000e+23 -4.460000e+00 8.200000e+01," &
        " 2.850000e+28 -5.470000e+00 3.003000e+03, 5.450000e+33" &
        " -6.610000e+00 7.559000e+03, 2.330000e+33 -6.110000e+00" &
        " 9.893000e+03, 3.270000e+25 -3.510000e+00 8.145000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_169 {
  chem_eq = "H + C4H8_1 --> H + C4H8_2"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 2.980000e+07 1.860000e+00" &
        " 3.575000e+03, 6.110000e+07 1.770000e+00 3.794000e+03," &
        " 4.780000e+09 1.240000e+00 5.152000e+03, 1.020000e+15" &
        " -2.500000e-01 9.233000e+03, 6.510000e+20 -1.820000e+00" &
        " 1.480600e+04, 4.440000e+19 -1.370000e+00 1.740900e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_169_rev {
  chem_eq = " H + C4H8_2 --> H + C4H8_1 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 2.980000e+07 1.860000e+00" &
        " 3.575000e+03, 6.110000e+07 1.770000e+00 3.794000e+03," &
        " 4.780000e+09 1.240000e+00 5.152000e+03, 1.020000e+15" &
        " -2.500000e-01 9.233000e+03, 6.510000e+20 -1.820000e+00" &
        " 1.480600e+04, 4.440000e+19 -1.370000e+00 1.740900e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_170 {
  chem_eq = "H + C4H8_1 --> H + C4H8_2"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.550000e+04 2.320000e+00" &
        " 7.049000e+03, 2.360000e+04 2.270000e+00 7.177000e+03," &
        " 6.600000e+05 1.860000e+00 8.201000e+03, 1.150000e+12" &
        " 1.100000e-01 1.278900e+04, 8.800000e+23 -3.170000e+00" &
        " 2.254600e+04, 3.720000e+31 -5.160000e+00 3.223400e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_170_rev {
  chem_eq = " H + C4H8_2 --> H + C4H8_1 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.550000e+04 2.320000e+00" &
        " 7.049000e+03, 2.360000e+04 2.270000e+00 7.177000e+03," &
        " 6.600000e+05 1.860000e+00 8.201000e+03, 1.150000e+12" &
        " 1.100000e-01 1.278900e+04, 8.800000e+23 -3.170000e+00" &
        " 2.254600e+04, 3.720000e+31 -5.160000e+00 3.223400e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_171 {
  chem_eq = "SC4H9 --> PC4H9"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 9.600000e+37 -1.104000e+01" &
        " 3.884000e+04, 6.050000e+40 -1.126000e+01 3.946100e+04," &
        " 1.640000e+47 -1.249000e+01 4.311200e+04, 6.530000e+55" &
        " -1.427000e+01 5.035100e+04, 2.130000e+56 -1.371000e+01" &
        " 5.486600e+04, 6.020000e+45 -1.007000e+01 5.339900e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_171_rev {
  chem_eq = " PC4H9 --> SC4H9 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 9.600000e+37 -1.104000e+01" &
        " 3.884000e+04, 6.050000e+40 -1.126000e+01 3.946100e+04," &
        " 1.640000e+47 -1.249000e+01 4.311200e+04, 6.530000e+55" &
        " -1.427000e+01 5.035100e+04, 2.130000e+56 -1.371000e+01" &
        " 5.486600e+04, 6.020000e+45 -1.007000e+01 5.339900e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_172 {
  chem_eq = "PC4H9 --> C2H5 + C2H4"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.440000e+34 -8.100000e+00" &
        " 2.839700e+04, 1.110000e+39 -9.050000e+00 3.189100e+04," &
        " 7.740000e+42 -9.780000e+00 3.577100e+04, 7.470000e+43" &
        " -9.670000e+00 3.872200e+04, 2.060000e+39 -7.970000e+00" &
        " 3.895500e+04, 1.480000e+29 -4.710000e+00 3.595000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_172_rev {
  chem_eq = " C2H5 + C2H4 --> PC4H9 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.440000e+34 -8.100000e+00" &
        " 2.839700e+04, 1.110000e+39 -9.050000e+00 3.189100e+04," &
        " 7.740000e+42 -9.780000e+00 3.577100e+04, 7.470000e+43" &
        " -9.670000e+00 3.872200e+04, 2.060000e+39 -7.970000e+00" &
        " 3.895500e+04, 1.480000e+29 -4.710000e+00 3.595000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_173 {
  chem_eq = "PC4H9 --> CH3 + C3H6"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.710000e+25 -5.810000e+00" &
        " 3.496500e+04, 1.850000e+27 -6.010000e+00 3.548100e+04," &
        " 2.460000e+32 -7.160000e+00 3.863700e+04, 2.050000e+42" &
        " -9.610000e+00 4.641500e+04, 4.980000e+48 -1.097000e+01" &
        " 5.445600e+04, 2.230000e+42 -8.680000e+00 5.660100e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_173_rev {
  chem_eq = " CH3 + C3H6 --> PC4H9 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 3.710000e+25 -5.810000e+00" &
        " 3.496500e+04, 1.850000e+27 -6.010000e+00 3.548100e+04," &
        " 2.460000e+32 -7.160000e+00 3.863700e+04, 2.050000e+42" &
        " -9.610000e+00 4.641500e+04, 4.980000e+48 -1.097000e+01" &
        " 5.445600e+04, 2.230000e+42 -8.680000e+00 5.660100e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_174 {
  chem_eq = "SC4H9 --> C2H5 + C2H4"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 8.300000e+25 -5.750000e+00" &
        " 3.934300e+04, 4.120000e+27 -5.940000e+00 3.985900e+04," &
        " 5.570000e+32 -7.100000e+00 4.302900e+04, 4.540000e+42" &
        " -9.540000e+00 5.083900e+04, 1.060000e+49 -1.090000e+01" &
        " 5.889900e+04, 9.940000e+42 -8.700000e+00 6.120300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_174_rev {
  chem_eq = " C2H5 + C2H4 --> SC4H9 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 8.300000e+25 -5.750000e+00" &
        " 3.934300e+04, 4.120000e+27 -5.940000e+00 3.985900e+04," &
        " 5.570000e+32 -7.100000e+00 4.302900e+04, 4.540000e+42" &
        " -9.540000e+00 5.083900e+04, 1.060000e+49 -1.090000e+01" &
        " 5.889900e+04, 9.940000e+42 -8.700000e+00 6.120300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_175 {
  chem_eq = "SC4H9 --> CH3 + C3H6"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 2.890000e+40 -9.760000e+00" &
        " 3.360100e+04, 1.800000e+44 -1.050000e+01 3.700700e+04," &
        " 2.510000e+46 -1.073000e+01 4.023700e+04, 4.740000e+44" &
        " -9.850000e+00 4.184100e+04, 3.790000e+37 -7.440000e+00" &
        " 4.060400e+04, 4.790000e+26 -4.010000e+00 3.689800e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_175_rev {
  chem_eq = " CH3 + C3H6 --> SC4H9 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-03" &
        " 1.000000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 2.890000e+40 -9.760000e+00" &
        " 3.360100e+04, 1.800000e+44 -1.050000e+01 3.700700e+04," &
        " 2.510000e+46 -1.073000e+01 4.023700e+04, 4.740000e+44" &
        " -9.850000e+00 4.184100e+04, 3.790000e+37 -7.440000e+00" &
        " 4.060400e+04, 4.790000e+26 -4.010000e+00 3.689800e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_176 {
  chem_eq = "C4H8_1 --> C4H8_2"
  arrhenius_coeff = 1.0000e+12 0.000 61000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_176_rev {
  chem_eq = " C4H8_2 --> C4H8_1 "
  arrhenius_coeff = 1.0000e+12 0.000 61000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_177 {
  chem_eq = "C4H8_1 --> H2 + C4H6"
  arrhenius_coeff = 3.0000e+13 0.000 70000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_177_rev {
  chem_eq = " H2 + C4H6 --> C4H8_1 "
  arrhenius_coeff = 3.0000e+13 0.000 70000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_178 {
  chem_eq = "C4H8_2 --> H2 + C4H6"
  arrhenius_coeff = 2.0000e+12 0.000 71000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_178_rev {
  chem_eq = " H2 + C4H6 --> C4H8_2 "
  arrhenius_coeff = 2.0000e+12 0.000 71000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_179 {
  chem_eq = "C4H6 --> H2 + C4H4"
  arrhenius_coeff = 2.5000e+15 0.000 94700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_179_rev {
  chem_eq = " H2 + C4H4 --> C4H6 "
  arrhenius_coeff = 2.5000e+15 0.000 94700.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_180 {
  chem_eq = "H + C4H6 --> C2H4 + C2H3"
  arrhenius_coeff = 6.4800e+32 -4.910 26478.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 8.640000e+33 -5.380000e+00 2.326400e+04, 6.480000e+32" &
        " -4.910000e+00 2.647800e+04, 1.170000e+20 -1.140000e+00" &
        " 2.302700e+04, 5.670000e+03 3.510000e+00 1.641500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_180_rev {
  chem_eq = " C2H4 + C2H3 --> H + C4H6 "
  arrhenius_coeff = 6.4800e+32 -4.910 26478.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 8.640000e+33 -5.380000e+00 2.326400e+04, 6.480000e+32" &
        " -4.910000e+00 2.647800e+04, 1.170000e+20 -1.140000e+00" &
        " 2.302700e+04, 5.670000e+03 3.510000e+00 1.641500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_181 {
  chem_eq = "H + C4H6 --> CH3 + C3H4_P"
  arrhenius_coeff = 2.0000e+12 0.000 7000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_181_rev {
  chem_eq = " CH3 + C3H4_P --> H + C4H6 "
  arrhenius_coeff = 2.0000e+12 0.000 7000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_182 {
  chem_eq = "H + C4H6 --> CH3 + C3H4_A"
  arrhenius_coeff = 2.0000e+12 0.000 7000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_182_rev {
  chem_eq = " CH3 + C3H4_A --> H + C4H6 "
  arrhenius_coeff = 2.0000e+12 0.000 7000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_183 {
  chem_eq = "C2H3 + C2H2 --> H + C4H4"
  arrhenius_coeff = 7.2000e+13 -0.480 6100.00
  press_rxn_param = PLOG "press_coeff: 1.320000e-02" &
        " 2.630000e-02 1.200000e-01 1.000000e+00 1.000000e+01" &
        " arrhenius_press: 7.200000e+13 -4.800000e-01 6.100000e+03," &
        " 5.000000e+14 -7.100000e-01 6.700000e+03, 4.600000e+16" &
        " -1.250000e+00 8.400000e+03, 2.000000e+18 -1.680000e+00" &
        " 1.060000e+04, 4.900000e+16 -1.130000e+00 1.180000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_183_rev {
  chem_eq = " H + C4H4 --> C2H3 + C2H2 "
  arrhenius_coeff = 7.2000e+13 -0.480 6100.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.320000e-02" &
        " 2.630000e-02 1.200000e-01 1.000000e+00 1.000000e+01" &
        " arrhenius_press: 7.200000e+13 -4.800000e-01 6.100000e+03," &
        " 5.000000e+14 -7.100000e-01 6.700000e+03, 4.600000e+16" &
        " -1.250000e+00 8.400000e+03, 2.000000e+18 -1.680000e+00" &
        " 1.060000e+04, 4.900000e+16 -1.130000e+00 1.180000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_184 {
  chem_eq = "2.0*C2H3 --> C4H6"
  arrhenius_coeff = 1.5000e+52 -11.970 16056.00
  press_rxn_param = PLOG "press_coeff: 2.630000e-02" &
        " 1.200000e-01 1.000000e+00 arrhenius_press: 7.000000e+57" &
        " -1.382000e+01 1.762900e+04, 1.500000e+52 -1.197000e+01" &
        " 1.605600e+04, 1.500000e+42 -8.840000e+00 1.248300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_184_rev {
  chem_eq = " C4H6 --> 2.0*C2H3 "
  arrhenius_coeff = 1.5000e+52 -11.970 16056.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 2.630000e-02" &
        " 1.200000e-01 1.000000e+00 arrhenius_press: 7.000000e+57" &
        " -1.382000e+01 1.762900e+04, 1.500000e+52 -1.197000e+01" &
        " 1.605600e+04, 1.500000e+42 -8.840000e+00 1.248300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_185 {
  chem_eq = "C2H2 + C2H --> H + C4H2"
  arrhenius_coeff = 9.6000e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_185_rev {
  chem_eq = " H + C4H2 --> C2H2 + C2H "
  arrhenius_coeff = 9.6000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_186 {
  chem_eq = "C4H6 --> H + C4H5"
  arrhenius_coeff = 5.7000e+36 -6.270 112353.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_187 {
  chem_eq = "C4H6 --> H + C4H5"
  arrhenius_coeff = 5.3000e+44 -8.620 123608.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_187_rev {
  chem_eq = " H + C4H5 --> C4H6 "
  arrhenius_coeff = 5.3000e+44 -8.620 123608.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_188 {
  chem_eq = "C2H3 + C4H6 --> C2H4 + C4H5"
  arrhenius_coeff = 5.0000e+13 0.000 22800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_188_rev {
  chem_eq = " C2H4 + C4H5 --> C2H3 + C4H6 "
  arrhenius_coeff = 5.0000e+13 0.000 22800.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_189 {
  chem_eq = "C3H3 + C4H6 --> C3H4_A + C4H5"
  arrhenius_coeff = 1.0000e+13 0.000 22500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_189_rev {
  chem_eq = " C3H4_A + C4H5 --> C3H3 + C4H6 "
  arrhenius_coeff = 1.0000e+13 0.000 22500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_190 {
  chem_eq = "C3H5_A + C4H6 --> C3H6 + C4H5"
  arrhenius_coeff = 1.0000e+13 0.000 22500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_190_rev {
  chem_eq = " C3H6 + C4H5 --> C3H5_A + C4H6 "
  arrhenius_coeff = 1.0000e+13 0.000 22500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_191 {
  chem_eq = "C2H3 + C2H2 --> C4H5"
  arrhenius_coeff = 1.1000e+31 -7.140 5600.00
  press_rxn_param = PLOG "press_coeff: 1.320000e-02" &
        " 2.630000e-02 1.200000e-01 1.000000e+00 1.000000e+01" &
        " arrhenius_press: 1.100000e+31 -7.140000e+00 5.600000e+03," &
        " 1.100000e+32 -7.330000e+00 6.200000e+03, 2.400000e+31" &
        " -6.950000e+00 5.600000e+03, 9.300000e+38 -8.760000e+00" &
        " 1.200000e+04, 8.100000e+37 -8.090000e+00 1.340000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_191_rev {
  chem_eq = " C4H5 --> C2H3 + C2H2 "
  arrhenius_coeff = 1.1000e+31 -7.140 5600.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.320000e-02" &
        " 2.630000e-02 1.200000e-01 1.000000e+00 1.000000e+01" &
        " arrhenius_press: 1.100000e+31 -7.140000e+00 5.600000e+03," &
        " 1.100000e+32 -7.330000e+00 6.200000e+03, 2.400000e+31" &
        " -6.950000e+00 5.600000e+03, 9.300000e+38 -8.760000e+00" &
        " 1.200000e+04, 8.100000e+37 -8.090000e+00 1.340000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_192 {
  chem_eq = "2.0*C2H3 --> H + C4H5"
  arrhenius_coeff = 1.1000e+24 -3.280 12395.00
  press_rxn_param = PLOG "press_coeff: 2.630000e-02" &
        " 1.200000e-01 1.000000e+00 arrhenius_press: 1.100000e+24" &
        " -3.280000e+00 1.239500e+04, 4.600000e+24 -3.380000e+00" &
        " 1.465000e+04, 2.400000e+20 -2.040000e+00 1.536100e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_192_rev {
  chem_eq = " H + C4H5 --> 2.0*C2H3 "
  arrhenius_coeff = 1.1000e+24 -3.280 12395.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 2.630000e-02" &
        " 1.200000e-01 1.000000e+00 arrhenius_press: 1.100000e+24" &
        " -3.280000e+00 1.239500e+04, 4.600000e+24 -3.380000e+00" &
        " 1.465000e+04, 2.400000e+20 -2.040000e+00 1.536100e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_193 {
  chem_eq = "H + C4H4 --> C4H5"
  arrhenius_coeff = 1.2000e+51 -12.570 12300.00
  press_rxn_param = PLOG "press_coeff: 1.320000e-02" &
        " 2.630000e-02 1.200000e-01 1.000000e+00 1.000000e+01" &
        " arrhenius_press: 1.200000e+51 -1.257000e+01 1.230000e+04," &
        " 4.200000e+50 -1.234000e+01 1.250000e+04, 1.100000e+50" &
        " -1.194000e+01 1.340000e+04, 1.300000e+51 -1.192000e+01" &
        " 1.650000e+04, 6.200000e+45 -1.008000e+01 1.580000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_193_rev {
  chem_eq = " C4H5 --> H + C4H4 "
  arrhenius_coeff = 1.2000e+51 -12.570 12300.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.320000e-02" &
        " 2.630000e-02 1.200000e-01 1.000000e+00 1.000000e+01" &
        " arrhenius_press: 1.200000e+51 -1.257000e+01 1.230000e+04," &
        " 4.200000e+50 -1.234000e+01 1.250000e+04, 1.100000e+50" &
        " -1.194000e+01 1.340000e+04, 1.300000e+51 -1.192000e+01" &
        " 1.650000e+04, 6.200000e+45 -1.008000e+01 1.580000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_194 {
  chem_eq = "H + C4H4 --> C4H5"
  arrhenius_coeff = 6.1000e+53 -13.190 14200.00
  press_rxn_param = PLOG "press_coeff: 1.320000e-02" &
        " 2.630000e-02 1.200000e-01 1.000000e+00 1.000000e+01" &
        " arrhenius_press: 6.100000e+53 -1.319000e+01 1.420000e+04," &
        " 9.600000e+52 -1.285000e+01 1.430000e+04, 2.100000e+52" &
        " -1.244000e+01 1.550000e+04, 4.900000e+51 -1.192000e+01" &
        " 1.770000e+04, 1.500000e+48 -1.058000e+01 1.880000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_194_rev {
  chem_eq = " C4H5 --> H + C4H4 "
  arrhenius_coeff = 6.1000e+53 -13.190 14200.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.320000e-02" &
        " 2.630000e-02 1.200000e-01 1.000000e+00 1.000000e+01" &
        " arrhenius_press: 6.100000e+53 -1.319000e+01 1.420000e+04," &
        " 9.600000e+52 -1.285000e+01 1.430000e+04, 2.100000e+52" &
        " -1.244000e+01 1.550000e+04, 4.900000e+51 -1.192000e+01" &
        " 1.770000e+04, 1.500000e+48 -1.058000e+01 1.880000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_195 {
  chem_eq = "H + C4H5 --> H2 + C4H4"
  arrhenius_coeff = 1.5000e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_195_rev {
  chem_eq = " H2 + C4H4 --> H + C4H5 "
  arrhenius_coeff = 1.5000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_196 {
  chem_eq = "C2H + C4H5 --> 2.0*C3H3"
  arrhenius_coeff = 4.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_196_rev {
  chem_eq = " 2.0*C3H3 --> C2H + C4H5 "
  arrhenius_coeff = 4.0000e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_197 {
  chem_eq = "H + C4H3 --> C4H4"
  arrhenius_coeff = 3.4000e+43 -9.010 12120.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_197_rev {
  chem_eq = " C4H4 --> H + C4H3 "
  arrhenius_coeff = 3.4000e+43 -9.010 12120.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_198 {
  chem_eq = "H + C4H4 --> H2 + C4H3"
  arrhenius_coeff = 6.6500e+05 2.530 12240.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_198_rev {
  chem_eq = " H2 + C4H3 --> H + C4H4 "
  arrhenius_coeff = 6.6500e+05 2.530 12240.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_199 {
  chem_eq = "CH3 + C4H4 --> CH4 + C4H3"
  arrhenius_coeff = 5.0000e+13 0.000 19800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_199_rev {
  chem_eq = " CH4 + C4H3 --> CH3 + C4H4 "
  arrhenius_coeff = 5.0000e+13 0.000 19800.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_200 {
  chem_eq = "CH3 + C4H4 --> C2H3 + C3H4_P"
  arrhenius_coeff = 2.7300e+09 1.080 19040.29
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 2.730000e+09" &
        " 1.080000e+00 1.904029e+04, 5.400000e+17 -1.180000e+00" &
        " 2.644381e+04, 3.660000e+17 -9.400000e-01 3.241374e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_201 {
  chem_eq = "C2H2 + C2H --> C4H3"
  arrhenius_coeff = 8.3000e+10 0.899 -363.00
  press_rxn_param = Troe_falloff "arrhenius_press: 1.24e+31" &
        " -4.718 1871.0 press_coeff: 1.000 100.0 5613. 1.339e+04"
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00" &
        " C2H2:2.50 C2H4:2.50 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_201_rev {
  chem_eq = " C4H3 --> C2H2 + C2H "
  arrhenius_coeff = 8.3000e+10 0.899 -363.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 1.24e+31" &
        " -4.718 1871.0 press_coeff: 1.000 100.0 5613. 1.339e+04"
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00" &
        " C2H2:2.50 C2H4:2.50 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_202 {
  chem_eq = "H + C4H2 --> C4H3"
  arrhenius_coeff = 1.1000e+42 -8.720 15300.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_202_rev {
  chem_eq = " C4H3 --> H + C4H2 "
  arrhenius_coeff = 1.1000e+42 -8.720 15300.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_203 {
  chem_eq = "H + C4H3 --> H2 + C4H2"
  arrhenius_coeff = 3.0000e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_203_rev {
  chem_eq = " H2 + C4H2 --> H + C4H3 "
  arrhenius_coeff = 3.0000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_204 {
  chem_eq = "C2H3 + C4H3 --> 2.0*C3H3"
  arrhenius_coeff = 4.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_204_rev {
  chem_eq = " 2.0*C3H3 --> C2H3 + C4H3 "
  arrhenius_coeff = 4.0000e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_205 {
  chem_eq = "2.0*C2H2 --> C2H3 + C2H"
  arrhenius_coeff = 1.7000e+14 0.000 92000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_205_rev {
  chem_eq = " C2H3 + C2H --> 2.0*C2H2 "
  arrhenius_coeff = 1.7000e+14 0.000 92000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_206 {
  chem_eq = "C2H6 + C2H2 --> C2H5 + C2H3"
  arrhenius_coeff = 2.0000e+14 0.000 60000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_206_rev {
  chem_eq = " C2H5 + C2H3 --> C2H6 + C2H2 "
  arrhenius_coeff = 2.0000e+14 0.000 60000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_207 {
  chem_eq = "2.0*CH3 --> H2 + C2H4"
  arrhenius_coeff = 5.0000e+14 0.000 32000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_208 {
  chem_eq = "C2H3 + C3H6 --> CH3 + C4H6"
  arrhenius_coeff = 1.8000e+11 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_209 {
  chem_eq = "C2H5 + IC3H7 --> C2H4 + C3H8"
  arrhenius_coeff = 5.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_210 {
  chem_eq = "CH3 + IC3H7 --> CH4 + C3H6"
  arrhenius_coeff = 6.0000e+10 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_211 {
  chem_eq = "C2H5 + C4H71_3 --> C2H6 + C4H6"
  arrhenius_coeff = 1.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_212 {
  chem_eq = "C2H2 + C3H4_A --> C2H + C3H5_A"
  arrhenius_coeff = 5.0000e+14 0.000 73000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_212_rev {
  chem_eq = " C2H + C3H5_A --> C2H2 + C3H4_A "
  arrhenius_coeff = 5.0000e+14 0.000 73000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_213 {
  chem_eq = "C2H2 + C3H6 --> C2H3 + C3H5_A"
  arrhenius_coeff = 4.0000e+13 0.000 46800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_213_rev {
  chem_eq = " C2H3 + C3H5_A --> C2H2 + C3H6 "
  arrhenius_coeff = 4.0000e+13 0.000 46800.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_214 {
  chem_eq = "C2H3 + IC4H7 --> C2H2 + IC4H8"
  arrhenius_coeff = 2.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_215 {
  chem_eq = "2.0*C3H3 --> C2H2 + C4H4"
  arrhenius_coeff = 1.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_216 {
  chem_eq = "2.0*C4H71_3 --> C4H8_2 + C4H6"
  arrhenius_coeff = 7.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_216_rev {
  chem_eq = " C4H8_2 + C4H6 --> 2.0*C4H71_3 "
  arrhenius_coeff = 7.0000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_217 {
  chem_eq = "IC4H7 + C4H71_3 --> IC4H8 + C4H6"
  arrhenius_coeff = 7.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_217_rev {
  chem_eq = " IC4H8 + C4H6 --> IC4H7 + C4H71_3 "
  arrhenius_coeff = 7.0000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_218 {
  chem_eq = "C2H3 + C4H71_3 --> C2H4 + C4H6"
  arrhenius_coeff = 2.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_218_rev {
  chem_eq = " C2H4 + C4H6 --> C2H3 + C4H71_3 "
  arrhenius_coeff = 2.0000e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_219 {
  chem_eq = "C4H5 --> C2H4 + C2H"
  arrhenius_coeff = 2.0000e+12 0.000 60000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_219_rev {
  chem_eq = " C2H4 + C2H --> C4H5 "
  arrhenius_coeff = 2.0000e+12 0.000 60000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_220 {
  chem_eq = "CH3 + C3H3 --> C4H6"
  arrhenius_coeff = 4.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_220_rev {
  chem_eq = " C4H6 --> CH3 + C3H3 "
  arrhenius_coeff = 4.0000e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_221 {
  chem_eq = "IC3H7 --> NC3H7"
  arrhenius_coeff = 2.0000e+12 0.000 42000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_221_rev {
  chem_eq = " NC3H7 --> IC3H7 "
  arrhenius_coeff = 2.0000e+12 0.000 42000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_222 {
  chem_eq = "CH3 + C4H6 --> C2H4 + C3H5_A"
  arrhenius_coeff = 1.0000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_223 {
  chem_eq = "C2H5 + C2H2 --> C4H71_4"
  arrhenius_coeff = 3.0000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_224 {
  chem_eq = "C2H5 + C3H6 --> C2H4 + NC3H7"
  arrhenius_coeff = 1.0000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_224_rev {
  chem_eq = " C2H4 + NC3H7 --> C2H5 + C3H6 "
  arrhenius_coeff = 1.0000e+11 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_225 {
  chem_eq = "C2H5 + C3H6 --> CH3 + C4H8_1"
  arrhenius_coeff = 9.0000e+10 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_225_rev {
  chem_eq = " CH3 + C4H8_1 --> C2H5 + C3H6 "
  arrhenius_coeff = 9.0000e+10 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_226 {
  chem_eq = "C2H5 + C4H6 --> C2H4 + C4H71_3"
  arrhenius_coeff = 1.0000e+10 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_227 {
  chem_eq = "C2H4 + IC3H7 --> C2H5 + C3H6"
  arrhenius_coeff = 1.0000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_227_rev {
  chem_eq = " C2H5 + C3H6 --> C2H4 + IC3H7 "
  arrhenius_coeff = 1.0000e+11 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_228 {
  chem_eq = "IC3H7 + C3H6 --> NC3H7 + C3H6"
  arrhenius_coeff = 5.0000e+10 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_228_rev {
  chem_eq = " NC3H7 + C3H6 --> IC3H7 + C3H6 "
  arrhenius_coeff = 5.0000e+10 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_229 {
  chem_eq = "C2H4 + C3H5_S --> CH3 + C4H6"
  arrhenius_coeff = 2.0000e+11 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_230 {
  chem_eq = "C3H6 + C3H5_S --> C2H5 + C4H6"
  arrhenius_coeff = 1.0000e+11 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_230_rev {
  chem_eq = " C2H5 + C4H6 --> C3H6 + C3H5_S "
  arrhenius_coeff = 1.0000e+11 0.000 6000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_231 {
  chem_eq = "C3H6 + C3H5_S --> C2H3 + C4H8_1"
  arrhenius_coeff = 7.0000e+10 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_231_rev {
  chem_eq = " C2H3 + C4H8_1 --> C3H6 + C3H5_S "
  arrhenius_coeff = 7.0000e+10 0.000 6000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_232 {
  chem_eq = "C2H4 + C4H71_3 --> 0.15*C2H5 + 0.85*C3H6 + 0.85*C3H5_A +" &
        " 0.15*C4H6"
  arrhenius_coeff = 2.0000e+11 0.000 13000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_233 {
  chem_eq = "C2H4 + C4H71_4 --> 0.8*C2H5 + 0.2*C3H6 + 0.1*C3H5_A +" &
        " 0.1*C3H5_S + 0.8*C4H6"
  arrhenius_coeff = 3.0000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_234 {
  chem_eq = "C3H6 + C4H71_4 --> NC3H7 + C4H6"
  arrhenius_coeff = 3.0000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_235 {
  chem_eq = "2.0*C4H6 --> C4H71_3 + C4H5"
  arrhenius_coeff = 5.0000e+14 0.000 64000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_235_rev {
  chem_eq = " C4H71_3 + C4H5 --> 2.0*C4H6 "
  arrhenius_coeff = 5.0000e+14 0.000 64000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_236 {
  chem_eq = "2.0*C3H6 --> NC3H7 + C3H5_A"
  arrhenius_coeff = 1.0000e+15 0.000 68000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_236_rev {
  chem_eq = " NC3H7 + C3H5_A --> 2.0*C3H6 "
  arrhenius_coeff = 1.0000e+15 0.000 68000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_237 {
  chem_eq = "2.0*C3H6 --> IC3H7 + C3H5_A"
  arrhenius_coeff = 1.0000e+15 0.000 68000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_237_rev {
  chem_eq = " IC3H7 + C3H5_A --> 2.0*C3H6 "
  arrhenius_coeff = 1.0000e+15 0.000 68000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_238 {
  chem_eq = "C2H6 --> H2 + C2H4"
  arrhenius_coeff = 3.0000e+13 0.000 78500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_238_rev {
  chem_eq = " H2 + C2H4 --> C2H6 "
  arrhenius_coeff = 3.0000e+13 0.000 78500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_239 {
  chem_eq = "C3H6 --> H2 + C3H4_A"
  arrhenius_coeff = 1.8000e+13 0.000 78000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_239_rev {
  chem_eq = " H2 + C3H4_A --> C3H6 "
  arrhenius_coeff = 1.8000e+13 0.000 78000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_240 {
  chem_eq = "C3H6 --> H2 + C3H4_P"
  arrhenius_coeff = 1.8000e+13 0.000 78000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_240_rev {
  chem_eq = " H2 + C3H4_P --> C3H6 "
  arrhenius_coeff = 1.8000e+13 0.000 78000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_241 {
  chem_eq = "C3H8 --> H2 + C3H6"
  arrhenius_coeff = 5.0000e+13 0.000 70000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_241_rev {
  chem_eq = " H2 + C3H6 --> C3H8 "
  arrhenius_coeff = 5.0000e+13 0.000 70000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_242 {
  chem_eq = "C4H4 --> H2 + C4H2"
  arrhenius_coeff = 5.0000e+11 0.000 66000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_242_rev {
  chem_eq = " H2 + C4H2 --> C4H4 "
  arrhenius_coeff = 5.0000e+11 0.000 66000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_243 {
  chem_eq = "C2H4 + C2H2 --> C4H6"
  arrhenius_coeff = 5.0000e+10 0.000 27000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_243_rev {
  chem_eq = " C4H6 --> C2H4 + C2H2 "
  arrhenius_coeff = 5.0000e+10 0.000 27000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_244 {
  chem_eq = "2.0*C2H4 --> C4H8_1"
  arrhenius_coeff = 1.0000e+10 0.000 40000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_244_rev {
  chem_eq = " C4H8_1 --> 2.0*C2H4 "
  arrhenius_coeff = 1.0000e+10 0.000 40000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_245 {
  chem_eq = "2.0*C2H2 --> C4H4"
  arrhenius_coeff = 5.0000e+11 0.000 37000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_245_rev {
  chem_eq = " C4H4 --> 2.0*C2H2 "
  arrhenius_coeff = 5.0000e+11 0.000 37000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_246 {
  chem_eq = "2.0*C2H2 --> H2 + C4H2"
  arrhenius_coeff = 1.0000e+16 0.000 68200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_246_rev {
  chem_eq = " H2 + C4H2 --> 2.0*C2H2 "
  arrhenius_coeff = 1.0000e+16 0.000 68200.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_247 {
  chem_eq = "2.0*C3H6 --> C2H4 + C4H8_2"
  arrhenius_coeff = 5.0000e+11 0.000 46000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_247_rev {
  chem_eq = " C2H4 + C4H8_2 --> 2.0*C3H6 "
  arrhenius_coeff = 5.0000e+11 0.000 46000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_248 {
  chem_eq = "H + C3H3 --> C3H4_A"
  arrhenius_coeff = 2.0500e+13 0.206 -173.00
  press_rxn_param = PLOG "press_coeff: 4.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 3.390000e+36 -7.410000e+00 6.337000e+03, 3.160000e+29" &
        " -5.000000e+00 4.711000e+03, 8.710000e+23 -3.200000e+00" &
        " 3.255000e+03, 2.050000e+13 2.060000e-01 -1.730000e+02"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_248_rev {
  chem_eq = " C3H4_A --> H + C3H3 "
  arrhenius_coeff = 2.0500e+13 0.206 -173.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 4.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 3.390000e+36 -7.410000e+00 6.337000e+03, 3.160000e+29" &
        " -5.000000e+00 4.711000e+03, 8.710000e+23 -3.200000e+00" &
        " 3.255000e+03, 2.050000e+13 2.060000e-01 -1.730000e+02"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_249 {
  chem_eq = "H + C3H3 --> C3H4_P"
  arrhenius_coeff = 6.4000e+13 0.102 -31.20
  press_rxn_param = PLOG "press_coeff: 4.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 3.630000e+36 -7.360000e+00 6.039000e+03, 7.940000e+29" &
        " -5.060000e+00 4.861000e+03, 1.070000e+24 -3.150000e+00" &
        " 3.261000e+03, 6.400000e+13 1.020000e-01 -3.120000e+01"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_249_rev {
  chem_eq = " C3H4_P --> H + C3H3 "
  arrhenius_coeff = 6.4000e+13 0.102 -31.20
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 4.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 3.630000e+36 -7.360000e+00 6.039000e+03, 7.940000e+29" &
        " -5.060000e+00 4.861000e+03, 1.070000e+24 -3.150000e+00" &
        " 3.261000e+03, 6.400000e+13 1.020000e-01 -3.120000e+01"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_250 {
  chem_eq = "H + C4H5 --> C2H4 + C2H2"
  arrhenius_coeff = 1.0000e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_250_rev {
  chem_eq = " C2H4 + C2H2 --> H + C4H5 "
  arrhenius_coeff = 1.0000e+13 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_251 {
  chem_eq = "C2H3 + C2H --> C4H4"
  arrhenius_coeff = 1.0000e+14 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_251_rev {
  chem_eq = " C4H4 --> C2H3 + C2H "
  arrhenius_coeff = 1.0000e+14 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_252 {
  chem_eq = "C2H2 + C4H4 --> C2H + C4H5"
  arrhenius_coeff = 1.0000e+14 0.000 95000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_252_rev {
  chem_eq = " C2H + C4H5 --> C2H2 + C4H4 "
  arrhenius_coeff = 1.0000e+14 0.000 95000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_253 {
  chem_eq = "C2H2 + C4H4 --> C2H3 + C4H3"
  arrhenius_coeff = 3.0000e+14 0.000 73000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_253_rev {
  chem_eq = " C2H3 + C4H3 --> C2H2 + C4H4 "
  arrhenius_coeff = 3.0000e+14 0.000 73000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_254 {
  chem_eq = "2.0*C4H4 --> C4H5 + C4H3"
  arrhenius_coeff = 5.0000e+16 0.000 81500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_254_rev {
  chem_eq = " C4H5 + C4H3 --> 2.0*C4H4 "
  arrhenius_coeff = 5.0000e+16 0.000 81500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_255 {
  chem_eq = "C3H3 + C4H71_3 --> C3H6 + C4H4"
  arrhenius_coeff = 4.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_256 {
  chem_eq = "C3H3 + IC4H7 --> C3H6 + C4H4"
  arrhenius_coeff = 4.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_257 {
  chem_eq = "2.0*C2H2 --> H + C4H3"
  arrhenius_coeff = 2.0000e+16 0.000 81500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_257_rev {
  chem_eq = " H + C4H3 --> 2.0*C2H2 "
  arrhenius_coeff = 2.0000e+16 0.000 81500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_258 {
  chem_eq = "C2H2 + C4H2 --> C2H + C4H3"
  arrhenius_coeff = 3.0000e+14 0.000 97000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_258_rev {
  chem_eq = " C2H + C4H3 --> C2H2 + C4H2 "
  arrhenius_coeff = 3.0000e+14 0.000 97000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_259 {
  chem_eq = "C2H + C3H4_P --> CH3 + C4H2"
  arrhenius_coeff = 3.0000e+11 0.000 8000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_260 {
  chem_eq = "2.0*C3H4_A --> C2H4 + C4H4"
  arrhenius_coeff = 5.0000e+11 0.000 29000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_261 {
  chem_eq = "2.0*C3H4_P --> C2H4 + C4H4"
  arrhenius_coeff = 5.0000e+11 0.000 29000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_262 {
  chem_eq = "H + C3H5_A --> C3H6"
  arrhenius_coeff = 1.0000e+14 0.000 0.00
  press_rxn_param = Troe_falloff "arrhenius_press: 1.33e+60" &
        " -12.000 5967.8 press_coeff: 0.02000 1097. 1.097e+04 6860."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_262_rev {
  chem_eq = " C3H6 --> H + C3H5_A "
  arrhenius_coeff = 1.0000e+14 0.000 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 1.33e+60" &
        " -12.000 5967.8 press_coeff: 0.02000 1097. 1.097e+04 6860."
  third_body_param = M_all "H2:2.00 CH4:2.00 C2H6:3.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_263 {
  chem_eq = "C3H6 + C3H5_A --> C2H4 + C4H71_4"
  arrhenius_coeff = 4.2000e+10 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_263_rev {
  chem_eq = " C2H4 + C4H71_4 --> C3H6 + C3H5_A "
  arrhenius_coeff = 4.2000e+10 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_264 {
  chem_eq = "C3H6 + C3H5_A --> H2 + CH3 + C5H6"
  arrhenius_coeff = 1.0000e+11 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_265 {
  chem_eq = "C3H6 + C3H5_A --> H2 + H + C5H5CH3"
  arrhenius_coeff = 4.0000e+10 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_266 {
  chem_eq = "C3H6 + C3H5_A --> C2H3 + C4H8_1"
  arrhenius_coeff = 2.4000e+10 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_266_rev {
  chem_eq = " C2H3 + C4H8_1 --> C3H6 + C3H5_A "
  arrhenius_coeff = 2.4000e+10 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_267 {
  chem_eq = "C3H5_A + C4H8_1 --> H2 + CH3 + C5H5CH3"
  arrhenius_coeff = 7.2000e+10 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_268 {
  chem_eq = "C3H5_A + C4H8_2 --> H2 + CH3 + C5H5CH3"
  arrhenius_coeff = 9.5000e+10 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_269 {
  chem_eq = "C3H5_A + C4H6 --> C2H5 + C5H6"
  arrhenius_coeff = 8.0000e+11 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_269_rev {
  chem_eq = " C2H5 + C5H6 --> C3H5_A + C4H6 "
  arrhenius_coeff = 8.0000e+11 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_270 {
  chem_eq = "C2H4 + PC4H9 --> NC3H7 + C3H6"
  arrhenius_coeff = 1.5000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_270_rev {
  chem_eq = " NC3H7 + C3H6 --> C2H4 + PC4H9 "
  arrhenius_coeff = 1.5000e+11 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_271 {
  chem_eq = "C3H5_A + C3H3 --> C2H2 + C4H6"
  arrhenius_coeff = 2.5000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_271_rev {
  chem_eq = " C2H2 + C4H6 --> C3H5_A + C3H3 "
  arrhenius_coeff = 2.5000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_272 {
  chem_eq = "C3H5_A + C4H71_3 --> C3H6 + C4H6"
  arrhenius_coeff = 5.0000e+11 0.000 3000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_272_rev {
  chem_eq = " C3H6 + C4H6 --> C3H5_A + C4H71_3 "
  arrhenius_coeff = 5.0000e+11 0.000 3000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_273 {
  chem_eq = "C3H5_A + IC4H7 --> C3H4_A + IC4H8"
  arrhenius_coeff = 2.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_273_rev {
  chem_eq = " C3H4_A + IC4H8 --> C3H5_A + IC4H7 "
  arrhenius_coeff = 2.0000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_274 {
  chem_eq = "C2H4 + IC4H7 --> 0.33*NC3H7 + 0.67*C3H6 + 0.67*C3H5_A +" &
        " 0.33*C3H4_A"
  arrhenius_coeff = 4.5000e+10 0.000 18000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_275 {
  chem_eq = "C3H4_P --> C3H4_A"
  arrhenius_coeff = 1.4000e+52 -10.860 95400.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 6.400000e+61 -1.459000e+01 8.820000e+04, 5.200000e+60" &
        " -1.393000e+01 9.110000e+04, 1.900000e+57 -1.262000e+01" &
        " 9.330000e+04, 1.400000e+52 -1.086000e+01 9.540000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_275_rev {
  chem_eq = " C3H4_A --> C3H4_P "
  arrhenius_coeff = 1.4000e+52 -10.860 95400.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 6.400000e+61 -1.459000e+01 8.820000e+04, 5.200000e+60" &
        " -1.393000e+01 9.110000e+04, 1.900000e+57 -1.262000e+01" &
        " 9.330000e+04, 1.400000e+52 -1.086000e+01 9.540000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_276 {
  chem_eq = "C2H3 + C4H2 --> C6H5"
  arrhenius_coeff = 5.0000e+11 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_277 {
  chem_eq = "C4H8_1 + C4H71_4 --> PC4H9 + C4H6"
  arrhenius_coeff = 3.0000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_277_rev {
  chem_eq = " PC4H9 + C4H6 --> C4H8_1 + C4H71_4 "
  arrhenius_coeff = 3.0000e+11 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_278 {
  chem_eq = "C4H8_2 + C4H71_4 --> SC4H9 + C4H6"
  arrhenius_coeff = 3.0000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_278_rev {
  chem_eq = " SC4H9 + C4H6 --> C4H8_2 + C4H71_4 "
  arrhenius_coeff = 3.0000e+11 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_279 {
  chem_eq = "NC4H10 --> H2 + 0.5*C4H8_1 + 0.5*C4H8_2"
  arrhenius_coeff = 5.0000e+13 0.000 71000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_280 {
  chem_eq = "IC4H10 --> H2 + IC4H8"
  arrhenius_coeff = 5.0000e+13 0.000 69000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_280_rev {
  chem_eq = " H2 + IC4H8 --> IC4H10 "
  arrhenius_coeff = 5.0000e+13 0.000 69000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_281 {
  chem_eq = "C2H2 + C3H4_A --> C5H6"
  arrhenius_coeff = 2.5000e+11 0.000 22000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_281_rev {
  chem_eq = " C5H6 --> C2H2 + C3H4_A "
  arrhenius_coeff = 2.5000e+11 0.000 22000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_282 {
  chem_eq = "C2H2 + C3H4_P --> C5H6"
  arrhenius_coeff = 2.5000e+10 0.000 22000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_282_rev {
  chem_eq = " C5H6 --> C2H2 + C3H4_P "
  arrhenius_coeff = 2.5000e+10 0.000 22000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_283 {
  chem_eq = "C2H3 + C4H6 --> CH3 + C5H6"
  arrhenius_coeff = 1.0000e+12 0.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_284 {
  chem_eq = "C2H4 + IC4H7 --> H + CH4 + C5H6"
  arrhenius_coeff = 1.0500e+11 0.000 18000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_285 {
  chem_eq = "C2H4 + C4H71_3 --> H2 + CH3 + C5H6"
  arrhenius_coeff = 1.5000e+11 0.000 13000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_286 {
  chem_eq = "CH3 + C4H3 --> C5H6"
  arrhenius_coeff = 1.0000e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_287 {
  chem_eq = "C5H5 --> C2H2 + C3H3"
  arrhenius_coeff = 6.8700e+55 -12.500 42065.01
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.780000e+78 -1.850000e+01 1.164720e+05," &
        " 1.130000e+74 -1.700000e+01 1.167270e+05, 1.640000e+66" &
        " -1.450000e+01 1.143410e+05, 1.350000e+57 -1.170000e+01" &
        " 1.115140e+05, 4.510000e+49 -9.400000e+00 1.120080e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_288 {
  chem_eq = "C2H2 + C3H3 --> C5H5"
  arrhenius_coeff = 6.8700e+55 -12.500 42065.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e-01 1.000000e+00 1.000000e+01 1.000000e+02" &
        " arrhenius_press: 7.510000e+68 -1.700000e+01 4.053600e+04," &
        " 1.050000e+64 -1.530000e+01 4.038000e+04, 8.480000e+54" &
        " -1.250000e+01 3.715300e+04, 1.870000e+45 -9.610000e+00" &
        " 3.396900e+04, 1.500000e+37 -7.100000e+00 3.406400e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_289 {
  chem_eq = "H + C5H5 --> C5H6"
  arrhenius_coeff = 2.5000e+13 0.280 -179.00
  press_rxn_param = Troe_falloff "arrhenius_press: 1.00e+81" &
        " -18.000 5000.0 press_coeff: 1.000 0.1000 585.0 6110."
  third_body_param = M_all "H2:2.00 CH4:2.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_289_rev {
  chem_eq = " C5H6 --> H + C5H5 "
  arrhenius_coeff = 2.5000e+13 0.280 -179.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 1.00e+81" &
        " -18.000 5000.0 press_coeff: 1.000 0.1000 585.0 6110."
  third_body_param = M_all "H2:2.00 CH4:2.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_290 {
  chem_eq = "CH3 + C5H5 --> C5H5CH3"
  arrhenius_coeff = 2.4600e+105 -27.028 47902.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 1.300000e+110 -2.918200e+01 4.190000e+04, 3.690000e+105" &
        " -2.702800e+01 4.790200e+04, 2.480000e+94 -2.347200e+01" &
        " 4.509300e+04, 8.490000e+69 -1.625400e+01 3.212700e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_290_rev {
  chem_eq = " C5H5CH3 --> CH3 + C5H5 "
  arrhenius_coeff = 2.4600e+105 -27.028 47902.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 1.300000e+110 -2.918200e+01 4.190000e+04, 3.690000e+105" &
        " -2.702800e+01 4.790200e+04, 2.480000e+94 -2.347200e+01" &
        " 4.509300e+04, 8.490000e+69 -1.625400e+01 3.212700e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_291 {
  chem_eq = "CH3 + C5H5 --> H2 + FULVENE"
  arrhenius_coeff = 9.2500e+01 2.190 8815.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 9.250000e+01 2.190000e+00 8.815000e+03, 1.140000e+15" &
        " -1.400000e+00 2.026600e+04, 3.620000e+31 -5.900000e+00" &
        " 3.563100e+04, 3.380000e+43 -9.050000e+00 4.974200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_292 {
  chem_eq = "H + C5H6 --> C2H2 + C3H5_A"
  arrhenius_coeff = 6.6500e+28 -3.910 29869.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 6.650000e+28" &
        " -3.910000e+00 2.986900e+04, 3.700000e+31 -4.600000e+00" &
        " 3.402500e+04, 6.860000e+32 -4.740000e+00 4.434700e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_292_rev {
  chem_eq = " C2H2 + C3H5_A --> H + C5H6 "
  arrhenius_coeff = 6.6500e+28 -3.910 29869.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 6.650000e+28" &
        " -3.910000e+00 2.986900e+04, 3.700000e+31 -4.600000e+00" &
        " 3.402500e+04, 6.860000e+32 -4.740000e+00 4.434700e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_293 {
  chem_eq = "H + C5H6 --> C2H4 + C3H3"
  arrhenius_coeff = 9.2300e+13 0.510 34563.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 9.230000e+13" &
        " 5.100000e-01 3.456300e+04, 5.320000e+23 -2.100000e+00" &
        " 4.489700e+04, 2.070000e+03 3.910000e+00 4.449200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_293_rev {
  chem_eq = " C2H4 + C3H3 --> H + C5H6 "
  arrhenius_coeff = 9.2300e+13 0.510 34563.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 9.230000e+13" &
        " 5.100000e-01 3.456300e+04, 5.320000e+23 -2.100000e+00" &
        " 4.489700e+04, 2.070000e+03 3.910000e+00 4.449200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_294 {
  chem_eq = "H + C5H6 --> H2 + C5H5"
  arrhenius_coeff = 5.7100e+05 2.390 2176.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_294_rev {
  chem_eq = " H2 + C5H5 --> H + C5H6 "
  arrhenius_coeff = 5.7100e+05 2.390 2176.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_295 {
  chem_eq = "CH3 + C5H6 --> CH4 + C5H5"
  arrhenius_coeff = 4.4453e+02 2.910 4727.71
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_295_rev {
  chem_eq = " CH4 + C5H5 --> CH3 + C5H6 "
  arrhenius_coeff = 4.4453e+02 2.910 4727.71
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_296 {
  chem_eq = "C2H3 + C5H6 --> C2H4 + C5H5"
  arrhenius_coeff = 1.4000e+06 2.000 4871.29
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_296_rev {
  chem_eq = " C2H4 + C5H5 --> C2H3 + C5H6 "
  arrhenius_coeff = 1.4000e+06 2.000 4871.29
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_297 {
  chem_eq = "C2H4 + C5H6 --> H + CH3 + C6H6"
  arrhenius_coeff = 7.5000e+10 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_298 {
  chem_eq = "C3H6 + C5H6 --> 2.0*CH3 + C6H6"
  arrhenius_coeff = 7.5000e+10 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_299 {
  chem_eq = "C4H8_1 + C5H6 --> CH3 + C2H5 + C6H6"
  arrhenius_coeff = 7.5000e+10 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_300 {
  chem_eq = "C4H8_2 + C5H6 --> CH3 + C2H5 + C6H6"
  arrhenius_coeff = 7.5000e+10 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_301 {
  chem_eq = "C4H6 + C5H6 --> CH3 + C2H3 + C6H6"
  arrhenius_coeff = 1.5000e+11 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_302 {
  chem_eq = "C4H5 + C5H6 --> C3H5_A + C6H6"
  arrhenius_coeff = 1.0000e+12 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_303 {
  chem_eq = "C5H6 + C5H5 --> C4H5 + C6H6"
  arrhenius_coeff = 2.0000e+13 0.000 25500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_304 {
  chem_eq = "C4H5 + C5H5 --> C3H4_A + C6H6"
  arrhenius_coeff = 2.0000e+12 0.000 3000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_305 {
  chem_eq = "C2H3 + C4H3 --> C6H6"
  arrhenius_coeff = 1.5000e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_306 {
  chem_eq = "C2H4 + C4H3 --> H + C6H6"
  arrhenius_coeff = 5.0000e+11 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_307 {
  chem_eq = "2.0*C4H4 --> C2H2 + C6H6"
  arrhenius_coeff = 2.5000e+14 0.000 44000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_308 {
  chem_eq = "C2H3 + C4H4 --> H + C6H6"
  arrhenius_coeff = 1.0000e+12 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_309 {
  chem_eq = "C2H + C4H4 --> C6H5"
  arrhenius_coeff = 3.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_309_rev {
  chem_eq = " C6H5 --> C2H + C4H4 "
  arrhenius_coeff = 3.0000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_310 {
  chem_eq = "C2H3 + C4H5 --> H2 + C6H6"
  arrhenius_coeff = 1.8400e-14 7.070 -3611.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_310_rev {
  chem_eq = " H2 + C6H6 --> C2H3 + C4H5 "
  arrhenius_coeff = 1.8400e-14 7.070 -3611.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_311 {
  chem_eq = "C3H4_A + C4H5 --> CH3 + C6H6"
  arrhenius_coeff = 2.5000e+11 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_312 {
  chem_eq = "C3H4_P + C4H5 --> CH3 + C6H6"
  arrhenius_coeff = 2.5000e+11 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_313 {
  chem_eq = "C4H5 + C4H2 --> C2H + C6H6"
  arrhenius_coeff = 5.0000e+11 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_314 {
  chem_eq = "C4H5 + C4H4 --> C2H3 + C6H6"
  arrhenius_coeff = 5.0000e+11 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_315 {
  chem_eq = "C2H3 + C4H6 --> H2 + H + C6H6"
  arrhenius_coeff = 5.6200e+11 0.000 3240.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_316 {
  chem_eq = "C2H4 + C4H6 --> 2.0*H2 + C6H6"
  arrhenius_coeff = 1.2500e+12 0.000 31000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_317 {
  chem_eq = "C2H + C4H5 --> C6H6"
  arrhenius_coeff = 1.0000e+13 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_318 {
  chem_eq = "C2H2 + C4H5 --> H + FULVENE"
  arrhenius_coeff = 1.2300e+20 -2.000 16152.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 2.500000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.520000e+15 -7.600000e-01" &
        " 8.767000e+03, 1.510000e-68 -7.600000e-01 8.767000e+03," &
        " 1.520000e+15 -7.600000e-01 8.769000e+03, 4.620000e+15" &
        " -8.900000e-01 9.142000e+03, 1.740000e+19 -1.860000e+00" &
        " 1.238300e+04, 1.230000e+20 -2.000000e+00 1.615200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_318_rev {
  chem_eq = " H + FULVENE --> C2H2 + C4H5 "
  arrhenius_coeff = 1.2300e+20 -2.000 16152.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 2.500000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.520000e+15 -7.600000e-01" &
        " 8.767000e+03, 1.510000e-68 -7.600000e-01 8.767000e+03," &
        " 1.520000e+15 -7.600000e-01 8.769000e+03, 4.620000e+15" &
        " -8.900000e-01 9.142000e+03, 1.740000e+19 -1.860000e+00" &
        " 1.238300e+04, 1.230000e+20 -2.000000e+00 1.615200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_319 {
  chem_eq = "C2H2 + C4H5 --> H + C6H6"
  arrhenius_coeff = 1.3700e+16 -1.000 8896.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 2.500000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.370000e+16 -1.000000e+00" &
        " 8.896000e+03, 2.940000e+16 -1.090000e+00 9.259000e+03," &
        " 1.370000e+16 -1.000000e+00 8.898000e+03, 1.380000e+16" &
        " -1.000000e+00 8.900000e+03, 1.690000e+16 -1.030000e+00" &
        " 8.967000e+03, 1.650000e+16 -1.010000e+00 9.480000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_319_rev {
  chem_eq = " H + C6H6 --> C2H2 + C4H5 "
  arrhenius_coeff = 1.3700e+16 -1.000 8896.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 2.500000e-02 1.000000e-01 1.000000e+00 1.000000e+01" &
        " 1.000000e+02 arrhenius_press: 1.370000e+16 -1.000000e+00" &
        " 8.896000e+03, 2.940000e+16 -1.090000e+00 9.259000e+03," &
        " 1.370000e+16 -1.000000e+00 8.898000e+03, 1.380000e+16" &
        " -1.000000e+00 8.900000e+03, 1.690000e+16 -1.030000e+00" &
        " 8.967000e+03, 1.650000e+16 -1.010000e+00 9.480000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_320 {
  chem_eq = "C2H2 + C4H3 --> C6H5"
  arrhenius_coeff = 9.6000e+70 -17.770 31300.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_320_rev {
  chem_eq = " C6H5 --> C2H2 + C4H3 "
  arrhenius_coeff = 9.6000e+70 -17.770 31300.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_321 {
  chem_eq = "C2H + C4H6 --> H + C6H6"
  arrhenius_coeff = 3.0000e+11 0.000 8000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_322 {
  chem_eq = "C2H2 + C4H4 --> C6H6"
  arrhenius_coeff = 4.5000e+11 0.000 30010.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_322_rev {
  chem_eq = " C6H6 --> C2H2 + C4H4 "
  arrhenius_coeff = 4.5000e+11 0.000 30010.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_323 {
  chem_eq = "2.0*C4H5 --> C2H4 + C6H6"
  arrhenius_coeff = 2.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_323_rev {
  chem_eq = " C2H4 + C6H6 --> 2.0*C4H5 "
  arrhenius_coeff = 2.0000e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_324 {
  chem_eq = "C3H4_A + C4H71_3 --> H2 + CH3 + C6H6"
  arrhenius_coeff = 2.0000e+12 0.000 18000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_325 {
  chem_eq = "C3H4_A + IC4H7 --> H2 + CH3 + C6H6"
  arrhenius_coeff = 2.0000e+12 0.000 18000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_326 {
  chem_eq = "C3H4_P + C4H71_3 --> H2 + CH3 + C6H6"
  arrhenius_coeff = 2.0000e+12 0.000 18000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_327 {
  chem_eq = "C3H4_P + IC4H7 --> H2 + CH3 + C6H6"
  arrhenius_coeff = 2.0000e+12 0.000 18000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_328 {
  chem_eq = "2.0*C3H3 --> H + C6H5"
  arrhenius_coeff = 1.0500e+54 -11.880 28757.00
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 7.900000e+53" &
        " -1.190000e+01 2.900000e+04, 8.500000e+47 -9.980000e+00" &
        " 3.680000e+04, 1.840000e+26 -3.880000e+00 2.900000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_328_rev {
  chem_eq = " H + C6H5 --> 2.0*C3H3 "
  arrhenius_coeff = 1.0500e+54 -11.880 28757.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 7.900000e+53" &
        " -1.190000e+01 2.900000e+04, 8.500000e+47 -9.980000e+00" &
        " 3.680000e+04, 1.840000e+26 -3.880000e+00 2.900000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_329 {
  chem_eq = "2.0*C3H3 --> C6H6"
  arrhenius_coeff = 6.4600e+26 -11.010 20320.00
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 1.290000e+69" &
        " -1.670000e+01 2.790000e+04, 3.160000e+55 -1.260000e+01" &
        " 2.230000e+04, 3.890000e+50 -1.100000e+01 2.030000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_329_rev {
  chem_eq = " C6H6 --> 2.0*C3H3 "
  arrhenius_coeff = 6.4600e+26 -11.010 20320.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 1.290000e+69" &
        " -1.670000e+01 2.790000e+04, 3.160000e+55 -1.260000e+01" &
        " 2.230000e+04, 3.890000e+50 -1.100000e+01 2.030000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_330 {
  chem_eq = "2.0*C3H3 --> FULVENE"
  arrhenius_coeff = 6.3100e+76 -19.070 31542.00
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 3.950000e-02 1.000000e+00 1.000000e+00 1.000000e+01" &
        " 1.000000e+01 arrhenius_press: 2.340000e+69 -1.700000e+01" &
        " 2.590000e+04, 1.730000e+44 -1.030000e+01 7.990000e+03," &
        " 5.490000e+62 -1.470000e+01 2.560000e+04, 3.020000e+35" &
        " -7.370000e+00 5.960000e+03, 5.620000e+60 -1.390000e+01" &
        " 2.710000e+04, 7.950000e+29 -5.500000e+00 4.670000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_330_rev {
  chem_eq = " FULVENE --> 2.0*C3H3 "
  arrhenius_coeff = 6.3100e+76 -19.070 31542.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 3.950000e-02 1.000000e+00 1.000000e+00 1.000000e+01" &
        " 1.000000e+01 arrhenius_press: 2.340000e+69 -1.700000e+01" &
        " 2.590000e+04, 1.730000e+44 -1.030000e+01 7.990000e+03," &
        " 5.490000e+62 -1.470000e+01 2.560000e+04, 3.020000e+35" &
        " -7.370000e+00 5.960000e+03, 5.620000e+60 -1.390000e+01" &
        " 2.710000e+04, 7.950000e+29 -5.500000e+00 4.670000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_331 {
  chem_eq = "C3H5_A + C3H3 --> 2.0*H + FULVENE"
  arrhenius_coeff = 2.6700e+35 -7.160 6810.01
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 2.670000e+35 -7.160000e+00 6.810010e+03, 6.230000e+38" &
        " -7.770000e+00 1.144290e+04, 8.200000e+12 -1.000000e-02" &
        " -1.471350e+03, 7.390000e+12 -0.000000e+00 -1.506870e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_332 {
  chem_eq = "C3H5_A + C3H3 --> 2.0*H + FULVENE"
  arrhenius_coeff = 7.7300e+41 -8.660 21473.95
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 7.730000e+41 -8.660000e+00 2.147395e+04, 1.920000e+38" &
        " -7.530000e+00 2.390330e+04, 4.790000e+23 -3.370000e+00" &
        " 1.822732e+04, 4.320000e+23 -3.360000e+00 1.819180e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_333 {
  chem_eq = "C6H6 --> H + C6H5"
  arrhenius_coeff = 5.5000e+38 -6.178 132000.00
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 1.350000e+108" &
        " -2.581000e+01 1.817500e+05, 6.310000e+60 -1.240000e+01" &
        " 1.480700e+05, 5.500000e+38 -6.178000e+00 1.320000e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_333_rev {
  chem_eq = " H + C6H5 --> C6H6 "
  arrhenius_coeff = 5.5000e+38 -6.178 132000.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 1.350000e+108" &
        " -2.581000e+01 1.817500e+05, 6.310000e+60 -1.240000e+01" &
        " 1.480700e+05, 5.500000e+38 -6.178000e+00 1.320000e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_334 {
  chem_eq = "FULVENE --> C6H6"
  arrhenius_coeff = 5.6200e+81 -19.360 121500.00
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 5.620000e+81" &
        " -1.936000e+01 1.215000e+05, 1.440000e+45 -8.900000e+00" &
        " 9.699900e+04, 2.950000e+31 -4.970000e+00 8.846500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_334_rev {
  chem_eq = " C6H6 --> FULVENE "
  arrhenius_coeff = 5.6200e+81 -19.360 121500.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 5.620000e+81" &
        " -1.936000e+01 1.215000e+05, 1.440000e+45 -8.900000e+00" &
        " 9.699900e+04, 2.950000e+31 -4.970000e+00 8.846500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_335 {
  chem_eq = "FULVENE --> H + C6H5"
  arrhenius_coeff = 2.5700e+97 -23.160 153470.00
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 2.570000e+97" &
        " -2.316000e+01 1.534700e+05, 2.240000e+68 -1.465000e+01" &
        " 1.425700e+05, 8.510000e+24 -4.970000e+00 1.133300e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_335_rev {
  chem_eq = " H + C6H5 --> FULVENE "
  arrhenius_coeff = 2.5700e+97 -23.160 153470.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 2.570000e+97" &
        " -2.316000e+01 1.534700e+05, 2.240000e+68 -1.465000e+01" &
        " 1.425700e+05, 8.510000e+24 -4.970000e+00 1.133300e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_336 {
  chem_eq = "C6H5 --> H + C2H2 + C4H2"
  arrhenius_coeff = 4.3000e+12 0.620 77300.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_337 {
  chem_eq = "H + FULVENE --> H + C6H6"
  arrhenius_coeff = 7.2600e+31 -4.748 18390.00
  press_rxn_param = PLOG "press_coeff: 1.320000e-01" &
        " 1.000000e+00 1.320000e+00 1.320000e+01 arrhenius_press:" &
        " 8.640000e+26 -3.446000e+00 1.264000e+04, 7.260000e+31" &
        " -4.748000e+00 1.839000e+04, 3.680000e+32 -4.931000e+00" &
        " 1.934000e+04, 6.080000e+36 -5.968000e+00 2.748000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_337_rev {
  chem_eq = " H + C6H6 --> H + FULVENE "
  arrhenius_coeff = 7.2600e+31 -4.748 18390.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.320000e-01" &
        " 1.000000e+00 1.320000e+00 1.320000e+01 arrhenius_press:" &
        " 8.640000e+26 -3.446000e+00 1.264000e+04, 7.260000e+31" &
        " -4.748000e+00 1.839000e+04, 3.680000e+32 -4.931000e+00" &
        " 1.934000e+04, 6.080000e+36 -5.968000e+00 2.748000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_338 {
  chem_eq = "H + C6H6 --> H2 + C6H5"
  arrhenius_coeff = 3.5900e+08 1.890 16052.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_338_rev {
  chem_eq = " H2 + C6H5 --> H + C6H6 "
  arrhenius_coeff = 3.5900e+08 1.890 16052.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_339 {
  chem_eq = "H + FULVENE --> H2 + C6H5"
  arrhenius_coeff = 3.5900e+08 1.890 16052.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_340 {
  chem_eq = "CH3 + C6H6 --> CH4 + C6H5"
  arrhenius_coeff = 5.0500e+02 3.190 13738.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_340_rev {
  chem_eq = " CH4 + C6H5 --> CH3 + C6H6 "
  arrhenius_coeff = 5.0500e+02 3.190 13738.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_341 {
  chem_eq = "CH3 + FULVENE --> CH4 + C6H5"
  arrhenius_coeff = 5.0500e+02 3.190 13738.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_342 {
  chem_eq = "C2H4 + C6H5 --> C2H3 + C6H6"
  arrhenius_coeff = 9.4500e-03 4.470 4471.80
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_342_rev {
  chem_eq = " C2H3 + C6H6 --> C2H4 + C6H5 "
  arrhenius_coeff = 9.4500e-03 4.470 4471.80
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_343 {
  chem_eq = "C2H2 + IC4H7 --> H + C5H5CH3"
  arrhenius_coeff = 1.0000e+13 0.000 14100.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_344 {
  chem_eq = "C2H2 + C4H71_3 --> H + C5H5CH3"
  arrhenius_coeff = 1.0000e+13 0.000 14100.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_345 {
  chem_eq = "C3H4_P + C3H4_A --> C5H5CH3"
  arrhenius_coeff = 8.0000e+10 0.000 21000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_345_rev {
  chem_eq = " C5H5CH3 --> C3H4_P + C3H4_A "
  arrhenius_coeff = 8.0000e+10 0.000 21000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_346 {
  chem_eq = "C5H5CH3 --> 2.0*H + FULVENE"
  arrhenius_coeff = 7.3500e+96 -24.412 116183.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 8.290000e+79 -2.073100e+01 9.503300e+04, 7.350000e+96" &
        " -2.441200e+01 1.161830e+05, 1.680000e+97 -2.398700e+01" &
        " 1.222730e+05, 2.870000e+86 -2.050900e+01 1.202310e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_347 {
  chem_eq = "C5H5CH3 --> 2.0*H + FULVENE"
  arrhenius_coeff = 1.7000e+110 -29.741 135237.00
  press_rxn_param = PLOG "press_coeff: 1.000000e+00" &
        " 1.000000e+01 1.000000e+02 arrhenius_press: 1.700000e+110" &
        " -2.974100e+01 1.352370e+05, 2.000000e+103 -2.613700e+01" &
        " 1.306300e+05, 4.300000e+93 -2.286900e+01 1.303190e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_348 {
  chem_eq = "C5H5CH3 --> H2 + FULVENE"
  arrhenius_coeff = 2.3300e+94 -25.069 118634.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 8.110000e+76 -2.149900e+01 9.668900e+04, 7.750000e+93" &
        " -2.506900e+01 1.186340e+05, 3.400000e+96 -2.522800e+01" &
        " 1.268500e+05, 8.310000e+83 -2.116100e+01 1.242630e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_348_rev {
  chem_eq = " H2 + FULVENE --> C5H5CH3 "
  arrhenius_coeff = 2.3300e+94 -25.069 118634.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+01 1.000000e+02 arrhenius_press:" &
        " 8.110000e+76 -2.149900e+01 9.668900e+04, 7.750000e+93" &
        " -2.506900e+01 1.186340e+05, 3.400000e+96 -2.522800e+01" &
        " 1.268500e+05, 8.310000e+83 -2.116100e+01 1.242630e+05"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_349 {
  chem_eq = "H + C5H5CH3 --> CH3 + C5H6"
  arrhenius_coeff = 1.5000e+13 0.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_350 {
  chem_eq = "H + C5H5CH3 --> H2 + H + FULVENE"
  arrhenius_coeff = 8.0000e+04 2.670 7645.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_351 {
  chem_eq = "H + C5H5CH3 --> H2 + H + C6H6"
  arrhenius_coeff = 1.2000e+07 2.000 2663.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_352 {
  chem_eq = "H + C5H5CH3 --> H2 + H + FULVENE"
  arrhenius_coeff = 2.7000e+04 2.670 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_353 {
  chem_eq = "CH3 + C5H5CH3 --> H + CH4 + FULVENE"
  arrhenius_coeff = 2.5000e+00 3.494 8599.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_354 {
  chem_eq = "CH3 + C5H5CH3 --> H + CH4 + C6H6"
  arrhenius_coeff = 3.6000e+05 2.000 3778.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_355 {
  chem_eq = "CH3 + C5H5CH3 --> H + CH4 + FULVENE"
  arrhenius_coeff = 8.8000e-01 3.494 1504.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_356 {
  chem_eq = "C5H5 + C5H5CH3 --> H + C6H6 + C5H6"
  arrhenius_coeff = 4.0000e+05 2.000 16500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_357 {
  chem_eq = "C5H5 + C5H5CH3 --> H + FULVENE + C5H6"
  arrhenius_coeff = 1.0000e+05 2.000 16500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_358 {
  chem_eq = "2.0*C5H6 --> 2.0*H2 + C10H8"
  arrhenius_coeff = 2.0000e+11 0.000 35000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_359 {
  chem_eq = "C4H2 + C6H6 --> C10H8"
  arrhenius_coeff = 1.5000e+11 0.000 22500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_360 {
  chem_eq = "C4H6 + C6H6 --> H2 + 2.0*H + C10H8"
  arrhenius_coeff = 1.0000e+11 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_361 {
  chem_eq = "C6H6 + C5H6 --> H + CH3 + C10H8"
  arrhenius_coeff = 2.0000e+11 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_362 {
  chem_eq = "C5H6 + C5H5 --> H2 + H + C10H8"
  arrhenius_coeff = 3.0000e+12 0.000 23000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_363 {
  chem_eq = "C6H6 + C5H5 --> CH3 + C10H8"
  arrhenius_coeff = 3.0000e+12 0.000 23000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_364 {
  chem_eq = "C4H3 + C6H6 --> H + C10H8"
  arrhenius_coeff = 5.0000e+11 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_365 {
  chem_eq = "C4H71_4 + C6H6 --> 2.0*H2 + H + C10H8"
  arrhenius_coeff = 4.5000e+11 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_366 {
  chem_eq = "C4H6 + C6H5 --> H2 + H + C10H8"
  arrhenius_coeff = 5.0000e+11 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_367 {
  chem_eq = "C4H4 + C6H5 --> H + C10H8"
  arrhenius_coeff = 1.2600e+04 2.610 1430.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_368 {
  chem_eq = "H + C10H8 --> C4H4 + C6H5"
  arrhenius_coeff = 5.0000e+12 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_369 {
  chem_eq = "H2 + C4H5 --> H + C4H6"
  arrhenius_coeff = 3.7800e+05 2.000 6910.63
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_370 {
  chem_eq = "H2 + C4H71_4 --> H + C4H8_1"
  arrhenius_coeff = 1.8900e+05 2.000 10425.49
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_371 {
  chem_eq = "H2 + IC4H7 --> H + IC4H8"
  arrhenius_coeff = 3.7800e+05 2.000 18239.88
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_372 {
  chem_eq = "H2 + LC5H7 --> H + LC5H8"
  arrhenius_coeff = 3.7800e+05 2.000 19648.56
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_373 {
  chem_eq = "CH4 + C4H5 --> CH3 + C4H6"
  arrhenius_coeff = 1.8900e+05 2.000 9418.91
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_374 {
  chem_eq = "CH4 + C4H71_4 --> CH3 + C4H8_1"
  arrhenius_coeff = 9.4500e+04 2.000 12868.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_375 {
  chem_eq = "CH4 + IC4H7 --> CH3 + IC4H8"
  arrhenius_coeff = 1.8900e+05 2.000 21138.84
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_376 {
  chem_eq = "CH4 + LC5H7 --> CH3 + LC5H8"
  arrhenius_coeff = 1.8900e+05 2.000 22616.91
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_377 {
  chem_eq = "CH4 + C7H7 --> CH3 + C7H8"
  arrhenius_coeff = 9.4500e+04 2.000 23748.43
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_378 {
  chem_eq = "CH4 + CH3C6H4 --> CH3 + C7H8"
  arrhenius_coeff = 2.5900e+08 1.000 6528.60
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_379 {
  chem_eq = "CH3 + C2H2 --> CH4 + C2H"
  arrhenius_coeff = 6.0000e+04 2.000 11196.53
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_380 {
  chem_eq = "C2H3 + C2H2 --> C2H4 + C2H"
  arrhenius_coeff = 1.3560e+05 2.000 9832.26
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_381 {
  chem_eq = "C2H5 + C2H2 --> C2H6 + C2H"
  arrhenius_coeff = 4.0000e+04 2.000 13265.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_382 {
  chem_eq = "C2H2 + NC3H7 --> C2H + C3H8"
  arrhenius_coeff = 2.7000e+04 2.000 13265.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_383 {
  chem_eq = "C2H2 + IC3H7 --> C2H + C3H8"
  arrhenius_coeff = 2.7000e+04 2.000 17064.96
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_384 {
  chem_eq = "C2H2 + C3H5_S --> C2H + C3H6"
  arrhenius_coeff = 5.4000e+04 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_385 {
  chem_eq = "C2H2 + C3H5_T --> C2H + C3H6"
  arrhenius_coeff = 5.4000e+04 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_386 {
  chem_eq = "C2H2 + C3H5_A --> C2H + C3H6"
  arrhenius_coeff = 1.0780e+05 2.000 20432.90
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_387 {
  chem_eq = "C2H2 + C4H5 --> C2H + C4H6"
  arrhenius_coeff = 5.4000e+04 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_388 {
  chem_eq = "C2H2 + C4H71_4 --> C2H + C4H8_1"
  arrhenius_coeff = 2.7000e+04 2.000 13144.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_389 {
  chem_eq = "C2H2 + C4H71_3 --> C2H + C4H8_2"
  arrhenius_coeff = 6.8000e+04 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_390 {
  chem_eq = "C2H2 + IC4H7 --> C2H + IC4H8"
  arrhenius_coeff = 5.4000e+04 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_391 {
  chem_eq = "C2H2 + C4H3 --> C2H + C4H4"
  arrhenius_coeff = 8.1800e+04 2.000 11442.22
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_392 {
  chem_eq = "C2H2 + C6H5 --> C2H + C6H6"
  arrhenius_coeff = 7.4000e+07 1.000 6728.39
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_393 {
  chem_eq = "C2H2 + C5H5 --> C2H + C5H6"
  arrhenius_coeff = 6.8000e+04 2.000 21810.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_394 {
  chem_eq = "C2H2 + LC5H7 --> C2H + LC5H8"
  arrhenius_coeff = 5.4000e+04 2.000 22952.95
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_395 {
  chem_eq = "C2H2 + C7H7 --> C2H + C7H8"
  arrhenius_coeff = 2.7000e+04 2.000 24090.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_396 {
  chem_eq = "C2H2 + CH3C6H4 --> C2H + C7H8"
  arrhenius_coeff = 7.4000e+07 1.000 6728.39
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_397 {
  chem_eq = "C2H6 + C2H3 --> C2H5 + C2H4"
  arrhenius_coeff = 6.1020e+05 2.000 5704.86
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_398 {
  chem_eq = "C2H6 + C3H5_S --> C2H5 + C3H6"
  arrhenius_coeff = 2.4300e+05 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_399 {
  chem_eq = "C2H6 + C3H5_T --> C2H5 + C3H6"
  arrhenius_coeff = 2.4300e+05 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_400 {
  chem_eq = "C2H6 + C4H5 --> C2H5 + C4H6"
  arrhenius_coeff = 2.4300e+05 2.000 6362.83
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_401 {
  chem_eq = "C2H6 + C4H71_4 --> C2H5 + C4H8_1"
  arrhenius_coeff = 1.2150e+05 2.000 9181.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_402 {
  chem_eq = "C2H6 + C4H71_3 --> C2H5 + C4H8_2"
  arrhenius_coeff = 3.0600e+05 2.000 16763.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_403 {
  chem_eq = "C2H6 + IC4H7 --> C2H5 + IC4H8"
  arrhenius_coeff = 2.4300e+05 2.000 16763.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_404 {
  chem_eq = "C2H6 + C4H3 --> C2H5 + C4H4"
  arrhenius_coeff = 3.6810e+05 2.000 7149.69
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_405 {
  chem_eq = "C2H6 + C6H5 --> C2H5 + C6H6"
  arrhenius_coeff = 3.3300e+08 1.000 3330.09
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_406 {
  chem_eq = "C2H6 + C5H5 --> C2H5 + C5H6"
  arrhenius_coeff = 3.0600e+05 2.000 17080.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_407 {
  chem_eq = "C2H6 + C3H3 --> C2H5 + C3H4_A"
  arrhenius_coeff = 2.4300e+05 2.000 16545.66
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_408 {
  chem_eq = "C2H6 + LC5H7 --> C2H5 + LC5H8"
  arrhenius_coeff = 2.4300e+05 2.000 18136.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_409 {
  chem_eq = "C2H6 + C7H7 --> C2H5 + C7H8"
  arrhenius_coeff = 1.2150e+05 2.000 19190.53
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_410 {
  chem_eq = "C2H6 + CH3C6H4 --> C2H5 + C7H8"
  arrhenius_coeff = 3.3300e+08 1.000 3864.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_411 {
  chem_eq = "C2H5 + C2H4 --> C2H6 + C2H3"
  arrhenius_coeff = 1.6000e+05 2.000 12804.86
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_412 {
  chem_eq = "C2H4 + C3H5_S --> C2H3 + C3H6"
  arrhenius_coeff = 2.1600e+05 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_413 {
  chem_eq = "C2H4 + C3H5_T --> C2H3 + C3H6"
  arrhenius_coeff = 2.1600e+05 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_414 {
  chem_eq = "C2H4 + C3H5_A --> C2H3 + C3H6"
  arrhenius_coeff = 4.3120e+05 2.000 20432.90
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_415 {
  chem_eq = "C2H4 + C4H71_4 --> C2H3 + C4H8_1"
  arrhenius_coeff = 1.0800e+05 2.000 13144.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_416 {
  chem_eq = "C2H4 + C4H71_3 --> C2H3 + C4H8_2"
  arrhenius_coeff = 2.7200e+05 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_417 {
  chem_eq = "C2H4 + IC4H7 --> C2H3 + IC4H8"
  arrhenius_coeff = 2.1600e+05 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_418 {
  chem_eq = "C2H4 + C4H3 --> C2H3 + C4H4"
  arrhenius_coeff = 3.2720e+05 2.000 10323.82
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_419 {
  chem_eq = "C2H4 + C3H3 --> C2H3 + C3H4_A"
  arrhenius_coeff = 2.1600e+05 2.000 20642.98
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_420 {
  chem_eq = "C2H4 + LC5H7 --> C2H3 + LC5H8"
  arrhenius_coeff = 2.1600e+05 2.000 22952.95
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_421 {
  chem_eq = "C2H4 + C7H7 --> C2H3 + C7H8"
  arrhenius_coeff = 1.0800e+05 2.000 24090.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_422 {
  chem_eq = "C2H4 + CH3C6H4 --> C2H3 + C7H8"
  arrhenius_coeff = 2.9600e+08 1.000 6728.39
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_423 {
  chem_eq = "C3H8 + C3H5_S --> NC3H7 + C3H6"
  arrhenius_coeff = 2.1600e+05 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_424 {
  chem_eq = "C3H8 + C3H5_T --> NC3H7 + C3H6"
  arrhenius_coeff = 2.1600e+05 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_425 {
  chem_eq = "C3H8 + C4H5 --> NC3H7 + C4H6"
  arrhenius_coeff = 2.1600e+05 2.000 6362.83
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_426 {
  chem_eq = "C3H8 + C4H71_4 --> NC3H7 + C4H8_1"
  arrhenius_coeff = 1.0800e+05 2.000 9181.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_427 {
  chem_eq = "C3H8 + C4H71_3 --> NC3H7 + C4H8_2"
  arrhenius_coeff = 2.7200e+05 2.000 16763.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_428 {
  chem_eq = "C3H8 + IC4H7 --> NC3H7 + IC4H8"
  arrhenius_coeff = 2.1600e+05 2.000 16763.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_429 {
  chem_eq = "C3H8 + C4H3 --> NC3H7 + C4H4"
  arrhenius_coeff = 3.2720e+05 2.000 7149.69
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_430 {
  chem_eq = "C3H8 + C6H5 --> NC3H7 + C6H6"
  arrhenius_coeff = 2.9600e+08 1.000 3330.09
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_431 {
  chem_eq = "C3H8 + C5H5 --> NC3H7 + C5H6"
  arrhenius_coeff = 2.7200e+05 2.000 17080.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_432 {
  chem_eq = "C3H8 + C3H3 --> NC3H7 + C3H4_A"
  arrhenius_coeff = 2.1600e+05 2.000 16545.66
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_433 {
  chem_eq = "C3H8 + LC5H7 --> NC3H7 + LC5H8"
  arrhenius_coeff = 2.1600e+05 2.000 18136.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_434 {
  chem_eq = "C3H8 + C7H7 --> NC3H7 + C7H8"
  arrhenius_coeff = 1.0800e+05 2.000 19190.53
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_435 {
  chem_eq = "C3H8 + CH3C6H4 --> NC3H7 + C7H8"
  arrhenius_coeff = 2.9600e+08 1.000 3864.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_436 {
  chem_eq = "C3H8 + C3H5_S --> IC3H7 + C3H6"
  arrhenius_coeff = 5.4000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_437 {
  chem_eq = "C3H8 + C3H5_T --> IC3H7 + C3H6"
  arrhenius_coeff = 5.4000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_438 {
  chem_eq = "C3H8 + C4H5 --> IC3H7 + C4H6"
  arrhenius_coeff = 5.4000e+04 2.000 4564.96
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_439 {
  chem_eq = "C3H8 + C4H71_4 --> IC3H7 + C4H8_1"
  arrhenius_coeff = 2.7000e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_440 {
  chem_eq = "C3H8 + C4H71_3 --> IC3H7 + C4H8_2"
  arrhenius_coeff = 6.8000e+04 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_441 {
  chem_eq = "C3H8 + IC4H7 --> IC3H7 + IC4H8"
  arrhenius_coeff = 5.4000e+04 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_442 {
  chem_eq = "C3H8 + C4H3 --> IC3H7 + C4H4"
  arrhenius_coeff = 8.1800e+04 2.000 5262.22
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_443 {
  chem_eq = "C3H8 + C6H5 --> IC3H7 + C6H6"
  arrhenius_coeff = 7.4000e+07 1.000 2064.96
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_444 {
  chem_eq = "C3H8 + C5H5 --> IC3H7 + C5H6"
  arrhenius_coeff = 6.8000e+04 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_445 {
  chem_eq = "C3H8 + C3H3 --> IC3H7 + C3H4_A"
  arrhenius_coeff = 5.4000e+04 2.000 14064.96
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_446 {
  chem_eq = "C3H8 + LC5H7 --> IC3H7 + LC5H8"
  arrhenius_coeff = 5.4000e+04 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_447 {
  chem_eq = "C3H8 + C7H7 --> IC3H7 + C7H8"
  arrhenius_coeff = 2.7000e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_448 {
  chem_eq = "C3H8 + CH3C6H4 --> IC3H7 + C7H8"
  arrhenius_coeff = 7.4000e+07 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_449 {
  chem_eq = "C2H3 + C3H6 --> C2H4 + C3H5_A"
  arrhenius_coeff = 2.0340e+05 2.000 4432.90
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_450 {
  chem_eq = "C3H6 + C3H5_S --> C3H6 + C3H5_A"
  arrhenius_coeff = 8.1000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_451 {
  chem_eq = "C3H6 + C3H5_T --> C3H6 + C3H5_A"
  arrhenius_coeff = 8.1000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_452 {
  chem_eq = "C3H6 + C4H71_4 --> C3H5_A + C4H8_1"
  arrhenius_coeff = 4.0500e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_453 {
  chem_eq = "C3H6 + C4H71_3 --> C3H5_A + C4H8_2"
  arrhenius_coeff = 1.0200e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_454 {
  chem_eq = "C3H6 + C4H3 --> C3H5_A + C4H4"
  arrhenius_coeff = 1.2270e+05 2.000 5710.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_455 {
  chem_eq = "C3H6 + C6H5 --> C3H5_A + C6H6"
  arrhenius_coeff = 1.1100e+08 1.000 2532.90
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_456 {
  chem_eq = "C3H6 + C5H5 --> C3H5_A + C5H6"
  arrhenius_coeff = 1.0200e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_457 {
  chem_eq = "C3H6 + LC5H7 --> C3H5_A + LC5H8"
  arrhenius_coeff = 8.1000e+04 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_458 {
  chem_eq = "C3H6 + C7H7 --> C3H5_A + C7H8"
  arrhenius_coeff = 4.0500e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_459 {
  chem_eq = "C3H6 + CH3C6H4 --> C3H5_A + C7H8"
  arrhenius_coeff = 1.1100e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_460 {
  chem_eq = "C2H3 + C3H6 --> C2H4 + C3H5_S"
  arrhenius_coeff = 2.7120e+05 2.000 8811.51
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_461 {
  chem_eq = "C2H5 + C3H6 --> C2H6 + C3H5_S"
  arrhenius_coeff = 8.0000e+04 2.000 12862.83
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_462 {
  chem_eq = "NC3H7 + C3H6 --> C3H8 + C3H5_S"
  arrhenius_coeff = 5.4000e+04 2.000 12862.83
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_463 {
  chem_eq = "IC3H7 + C3H6 --> C3H8 + C3H5_S"
  arrhenius_coeff = 5.4000e+04 2.000 16662.41
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_464 {
  chem_eq = "C3H6 + C3H5_T --> C3H6 + C3H5_S"
  arrhenius_coeff = 1.0800e+05 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_465 {
  chem_eq = "C3H6 + C3H5_A --> C3H6 + C3H5_S"
  arrhenius_coeff = 2.1560e+05 2.000 20432.90
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_466 {
  chem_eq = "C3H6 + C4H5 --> C3H5_S + C4H6"
  arrhenius_coeff = 1.0800e+05 2.000 9558.50
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_467 {
  chem_eq = "C3H6 + C4H71_4 --> C3H5_S + C4H8_1"
  arrhenius_coeff = 5.4000e+04 2.000 13144.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_468 {
  chem_eq = "C3H6 + C4H71_3 --> C3H5_S + C4H8_2"
  arrhenius_coeff = 1.3600e+05 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_469 {
  chem_eq = "C3H6 + C4H3 --> C3H5_S + C4H4"
  arrhenius_coeff = 1.6360e+05 2.000 10464.60
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_470 {
  chem_eq = "C3H6 + C6H5 --> C3H5_S + C6H6"
  arrhenius_coeff = 1.4800e+08 1.000 5707.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_471 {
  chem_eq = "C3H6 + C5H5 --> C3H5_S + C5H6"
  arrhenius_coeff = 1.3600e+05 2.000 21810.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_472 {
  chem_eq = "C3H6 + C3H3 --> C3H5_S + C3H4_A"
  arrhenius_coeff = 1.0800e+05 2.000 20789.97
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_473 {
  chem_eq = "C3H6 + LC5H7 --> C3H5_S + LC5H8"
  arrhenius_coeff = 1.0800e+05 2.000 22952.95
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_474 {
  chem_eq = "C3H6 + C7H7 --> C3H5_S + C7H8"
  arrhenius_coeff = 5.4000e+04 2.000 24090.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_475 {
  chem_eq = "C3H6 + CH3C6H4 --> C3H5_S + C7H8"
  arrhenius_coeff = 1.4800e+08 1.000 6728.39
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_476 {
  chem_eq = "C2H3 + C3H6 --> C2H4 + C3H5_T"
  arrhenius_coeff = 9.4920e+04 2.000 8811.51
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_477 {
  chem_eq = "C2H5 + C3H6 --> C2H6 + C3H5_T"
  arrhenius_coeff = 2.8000e+04 2.000 12862.83
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_478 {
  chem_eq = "NC3H7 + C3H6 --> C3H8 + C3H5_T"
  arrhenius_coeff = 1.8900e+04 2.000 12862.83
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_479 {
  chem_eq = "IC3H7 + C3H6 --> C3H8 + C3H5_T"
  arrhenius_coeff = 1.8900e+04 2.000 16662.41
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_480 {
  chem_eq = "C3H6 + C3H5_S --> C3H6 + C3H5_T"
  arrhenius_coeff = 3.7800e+04 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_481 {
  chem_eq = "C3H6 + C3H5_A --> C3H6 + C3H5_T"
  arrhenius_coeff = 7.5460e+04 2.000 20432.90
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_482 {
  chem_eq = "C3H6 + C4H5 --> C3H5_T + C4H6"
  arrhenius_coeff = 3.7800e+04 2.000 9558.50
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_483 {
  chem_eq = "C3H6 + C4H71_4 --> C3H5_T + C4H8_1"
  arrhenius_coeff = 1.8900e+04 2.000 13144.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_484 {
  chem_eq = "C3H6 + C4H71_3 --> C3H5_T + C4H8_2"
  arrhenius_coeff = 4.7600e+04 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_485 {
  chem_eq = "C3H6 + C4H3 --> C3H5_T + C4H4"
  arrhenius_coeff = 5.7260e+04 2.000 10464.60
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_486 {
  chem_eq = "C3H6 + C6H5 --> C3H5_T + C6H6"
  arrhenius_coeff = 5.1800e+07 1.000 5707.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_487 {
  chem_eq = "C3H6 + C5H5 --> C3H5_T + C5H6"
  arrhenius_coeff = 4.7600e+04 2.000 21810.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_488 {
  chem_eq = "C3H6 + C3H3 --> C3H5_T + C3H4_A"
  arrhenius_coeff = 3.7800e+04 2.000 20789.97
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_489 {
  chem_eq = "C3H6 + LC5H7 --> C3H5_T + LC5H8"
  arrhenius_coeff = 3.7800e+04 2.000 22952.95
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_490 {
  chem_eq = "C3H6 + C7H7 --> C3H5_T + C7H8"
  arrhenius_coeff = 1.8900e+04 2.000 24090.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_491 {
  chem_eq = "C3H6 + CH3C6H4 --> C3H5_T + C7H8"
  arrhenius_coeff = 5.1800e+07 1.000 6728.39
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_492 {
  chem_eq = "C2H5 + C3H4_A --> C2H6 + C3H3"
  arrhenius_coeff = 1.2000e+06 2.000 13545.66
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_493 {
  chem_eq = "NC3H7 + C3H4_A --> C3H8 + C3H3"
  arrhenius_coeff = 8.1000e+05 2.000 13545.66
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_494 {
  chem_eq = "IC3H7 + C3H4_A --> C3H8 + C3H3"
  arrhenius_coeff = 8.1000e+05 2.000 17345.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_495 {
  chem_eq = "C3H5_S + C3H4_A --> C3H6 + C3H3"
  arrhenius_coeff = 1.6200e+06 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_496 {
  chem_eq = "C3H5_T + C3H4_A --> C3H6 + C3H3"
  arrhenius_coeff = 1.6200e+06 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_497 {
  chem_eq = "C3H4_A + C4H71_4 --> C3H3 + C4H8_1"
  arrhenius_coeff = 8.1000e+05 2.000 13144.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_498 {
  chem_eq = "C3H4_A + C4H71_3 --> C3H3 + C4H8_2"
  arrhenius_coeff = 2.0400e+06 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_499 {
  chem_eq = "C3H4_A + IC4H7 --> C3H3 + IC4H8"
  arrhenius_coeff = 1.6200e+06 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_500 {
  chem_eq = "C3H4_A + C4H3 --> C3H3 + C4H4"
  arrhenius_coeff = 2.4540e+06 2.000 12122.91
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_501 {
  chem_eq = "C3H4_A + C6H5 --> C3H3 + C6H6"
  arrhenius_coeff = 2.2200e+09 1.000 7439.11
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_502 {
  chem_eq = "C3H4_A + C5H5 --> C3H3 + C5H6"
  arrhenius_coeff = 2.0400e+06 2.000 21810.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_503 {
  chem_eq = "C3H4_A + LC5H7 --> C3H3 + LC5H8"
  arrhenius_coeff = 1.6200e+06 2.000 22952.95
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_504 {
  chem_eq = "C3H4_A + C7H7 --> C3H3 + C7H8"
  arrhenius_coeff = 8.1000e+05 2.000 24090.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_505 {
  chem_eq = "C3H4_A + CH3C6H4 --> C3H3 + C7H8"
  arrhenius_coeff = 2.2200e+09 1.000 6728.39
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_506 {
  chem_eq = "C2H5 + C3H4_P --> C2H6 + C3H3"
  arrhenius_coeff = 1.2000e+06 2.000 13545.66
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_507 {
  chem_eq = "NC3H7 + C3H4_P --> C3H8 + C3H3"
  arrhenius_coeff = 8.1000e+05 2.000 13545.66
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_508 {
  chem_eq = "IC3H7 + C3H4_P --> C3H8 + C3H3"
  arrhenius_coeff = 8.1000e+05 2.000 17345.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_509 {
  chem_eq = "C3H5_S + C3H4_P --> C3H6 + C3H3"
  arrhenius_coeff = 1.6200e+06 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_510 {
  chem_eq = "C3H5_T + C3H4_P --> C3H6 + C3H3"
  arrhenius_coeff = 1.6200e+06 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_511 {
  chem_eq = "C3H4_P + C4H5 --> C3H3 + C4H6"
  arrhenius_coeff = 1.6200e+06 2.000 11289.97
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_512 {
  chem_eq = "C3H4_P + C4H71_4 --> C3H3 + C4H8_1"
  arrhenius_coeff = 8.1000e+05 2.000 13144.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_513 {
  chem_eq = "C3H4_P + C4H71_3 --> C3H3 + C4H8_2"
  arrhenius_coeff = 2.0400e+06 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_514 {
  chem_eq = "C3H4_P + IC4H7 --> C3H3 + IC4H8"
  arrhenius_coeff = 1.6200e+06 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_515 {
  chem_eq = "C3H4_P + C4H3 --> C3H3 + C4H4"
  arrhenius_coeff = 2.4540e+06 2.000 12122.91
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_516 {
  chem_eq = "C3H4_P + C6H5 --> C3H3 + C6H6"
  arrhenius_coeff = 2.2200e+09 1.000 7439.11
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_517 {
  chem_eq = "C3H4_P + C5H5 --> C3H3 + C5H6"
  arrhenius_coeff = 2.0400e+06 2.000 21810.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_518 {
  chem_eq = "C3H4_P + LC5H7 --> C3H3 + LC5H8"
  arrhenius_coeff = 1.6200e+06 2.000 22952.95
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_519 {
  chem_eq = "C3H4_P + C7H7 --> C3H3 + C7H8"
  arrhenius_coeff = 8.1000e+05 2.000 24090.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_520 {
  chem_eq = "C3H4_P + CH3C6H4 --> C3H3 + C7H8"
  arrhenius_coeff = 2.2200e+09 1.000 6728.39
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_521 {
  chem_eq = "C2H3 + C4H4 --> C2H4 + C4H3"
  arrhenius_coeff = 9.4920e+05 2.000 8723.82
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_522 {
  chem_eq = "C2H5 + C4H4 --> C2H6 + C4H3"
  arrhenius_coeff = 2.8000e+05 2.000 12649.68
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_523 {
  chem_eq = "NC3H7 + C4H4 --> C3H8 + C4H3"
  arrhenius_coeff = 1.8900e+05 2.000 12649.68
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_524 {
  chem_eq = "IC3H7 + C4H4 --> C3H8 + C4H3"
  arrhenius_coeff = 1.8900e+05 2.000 16423.94
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_525 {
  chem_eq = "C3H5_S + C4H4 --> C3H6 + C4H3"
  arrhenius_coeff = 3.7800e+05 2.000 10322.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_526 {
  chem_eq = "C3H5_T + C4H4 --> C3H6 + C4H3"
  arrhenius_coeff = 3.7800e+05 2.000 10322.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_527 {
  chem_eq = "C3H5_A + C4H4 --> C3H6 + C4H3"
  arrhenius_coeff = 7.5460e+05 2.000 20110.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_528 {
  chem_eq = "C4H5 + C4H4 --> C4H6 + C4H3"
  arrhenius_coeff = 3.7800e+05 2.000 9464.60
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_529 {
  chem_eq = "C4H71_4 + C4H4 --> C4H8_1 + C4H3"
  arrhenius_coeff = 1.8900e+05 2.000 12868.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_530 {
  chem_eq = "C4H71_3 + C4H4 --> C4H8_2 + C4H3"
  arrhenius_coeff = 4.7600e+05 2.000 21138.84
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_531 {
  chem_eq = "IC4H7 + C4H4 --> IC4H8 + C4H3"
  arrhenius_coeff = 3.7800e+05 2.000 21138.84
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_532 {
  chem_eq = "C4H4 + C6H5 --> C4H3 + C6H6"
  arrhenius_coeff = 5.1800e+08 1.000 5670.82
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_533 {
  chem_eq = "C4H4 + C5H5 --> C4H3 + C5H6"
  arrhenius_coeff = 4.7600e+05 2.000 21480.69
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_534 {
  chem_eq = "C3H3 + C4H4 --> C3H4_A + C4H3"
  arrhenius_coeff = 3.7800e+05 2.000 20622.91
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_535 {
  chem_eq = "C4H4 + LC5H7 --> C4H3 + LC5H8"
  arrhenius_coeff = 3.7800e+05 2.000 22616.91
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_536 {
  chem_eq = "C4H4 + C7H7 --> C4H3 + C7H8"
  arrhenius_coeff = 1.8900e+05 2.000 23748.43
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_537 {
  chem_eq = "C4H4 + CH3C6H4 --> C4H3 + C7H8"
  arrhenius_coeff = 5.1800e+08 1.000 6528.60
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_538 {
  chem_eq = "H + C4H6 --> H2 + C4H5"
  arrhenius_coeff = 7.2000e+07 2.000 9360.63
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_539 {
  chem_eq = "CH3 + C4H6 --> CH4 + C4H5"
  arrhenius_coeff = 7.2000e+05 2.000 10218.91
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_540 {
  chem_eq = "C2H5 + C4H6 --> C2H6 + C4H5"
  arrhenius_coeff = 4.8000e+05 2.000 12862.83
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_541 {
  chem_eq = "NC3H7 + C4H6 --> C3H8 + C4H5"
  arrhenius_coeff = 3.2400e+05 2.000 12862.83
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_542 {
  chem_eq = "IC3H7 + C4H6 --> C3H8 + C4H5"
  arrhenius_coeff = 3.2400e+05 2.000 16662.41
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_543 {
  chem_eq = "C3H5_S + C4H6 --> C3H6 + C4H5"
  arrhenius_coeff = 6.4800e+05 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_544 {
  chem_eq = "C3H5_T + C4H6 --> C3H6 + C4H5"
  arrhenius_coeff = 6.4800e+05 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_545 {
  chem_eq = "C4H71_4 + C4H6 --> C4H8_1 + C4H5"
  arrhenius_coeff = 3.2400e+05 2.000 13144.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_546 {
  chem_eq = "C4H71_3 + C4H6 --> C4H8_2 + C4H5"
  arrhenius_coeff = 8.1600e+05 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_547 {
  chem_eq = "IC4H7 + C4H6 --> IC4H8 + C4H5"
  arrhenius_coeff = 6.4800e+05 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_548 {
  chem_eq = "C4H6 + C4H3 --> C4H5 + C4H4"
  arrhenius_coeff = 9.8160e+05 2.000 10464.60
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_549 {
  chem_eq = "C4H6 + C6H5 --> C4H5 + C6H6"
  arrhenius_coeff = 8.8800e+08 1.000 5707.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_550 {
  chem_eq = "C4H6 + C5H5 --> C4H5 + C5H6"
  arrhenius_coeff = 8.1600e+05 2.000 21810.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_551 {
  chem_eq = "C4H6 + LC5H7 --> C4H5 + LC5H8"
  arrhenius_coeff = 6.4800e+05 2.000 22952.95
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_552 {
  chem_eq = "C4H6 + C7H7 --> C4H5 + C7H8"
  arrhenius_coeff = 3.2400e+05 2.000 24090.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_553 {
  chem_eq = "C4H6 + CH3C6H4 --> C4H5 + C7H8"
  arrhenius_coeff = 8.8800e+08 1.000 6728.39
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_554 {
  chem_eq = "H + C4H8_1 --> H2 + C4H71_3"
  arrhenius_coeff = 1.8000e+07 2.000 4389.88
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_555 {
  chem_eq = "CH3 + C4H8_1 --> CH4 + C4H71_3"
  arrhenius_coeff = 1.8000e+05 2.000 5638.84
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_556 {
  chem_eq = "C2H3 + C4H8_1 --> C2H4 + C4H71_3"
  arrhenius_coeff = 4.0680e+05 2.000 4567.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_557 {
  chem_eq = "C2H5 + C4H8_1 --> C2H6 + C4H71_3"
  arrhenius_coeff = 1.2000e+05 2.000 6963.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_558 {
  chem_eq = "NC3H7 + C4H8_1 --> C3H8 + C4H71_3"
  arrhenius_coeff = 8.1000e+04 2.000 6963.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_559 {
  chem_eq = "IC3H7 + C4H8_1 --> C3H8 + C4H71_3"
  arrhenius_coeff = 8.1000e+04 2.000 10163.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_560 {
  chem_eq = "C3H5_S + C4H8_1 --> C3H6 + C4H71_3"
  arrhenius_coeff = 1.6200e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_561 {
  chem_eq = "C3H5_T + C4H8_1 --> C3H6 + C4H71_3"
  arrhenius_coeff = 1.6200e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_562 {
  chem_eq = "C3H5_A + C4H8_1 --> C3H6 + C4H71_3"
  arrhenius_coeff = 3.2340e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_563 {
  chem_eq = "C4H8_1 + C4H5 --> C4H71_3 + C4H6"
  arrhenius_coeff = 1.6200e+05 2.000 5167.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_564 {
  chem_eq = "C4H8_1 + C4H71_4 --> C4H8_1 + C4H71_3"
  arrhenius_coeff = 8.1000e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_565 {
  chem_eq = "C4H8_1 + C4H71_3 --> C4H8_2 + C4H71_3"
  arrhenius_coeff = 2.0400e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_566 {
  chem_eq = "IC4H7 + C4H8_1 --> IC4H8 + C4H71_3"
  arrhenius_coeff = 1.6200e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_567 {
  chem_eq = "C4H8_1 + C4H3 --> C4H71_3 + C4H4"
  arrhenius_coeff = 2.4540e+05 2.000 5838.84
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_568 {
  chem_eq = "C4H8_1 + C6H5 --> C4H71_3 + C6H6"
  arrhenius_coeff = 2.2200e+08 1.000 2667.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_569 {
  chem_eq = "C4H8_1 + C5H5 --> C4H71_3 + C5H6"
  arrhenius_coeff = 2.0400e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_570 {
  chem_eq = "C3H3 + C4H8_1 --> C3H4_A + C4H71_3"
  arrhenius_coeff = 1.6200e+05 2.000 14667.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_571 {
  chem_eq = "C4H8_1 + LC5H7 --> C4H71_3 + LC5H8"
  arrhenius_coeff = 1.6200e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_572 {
  chem_eq = "C4H8_1 + C7H7 --> C4H71_3 + C7H8"
  arrhenius_coeff = 8.1000e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_573 {
  chem_eq = "C4H8_1 + CH3C6H4 --> C4H71_3 + C7H8"
  arrhenius_coeff = 2.2200e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_574 {
  chem_eq = "H + C4H8_1 --> H2 + C4H71_4"
  arrhenius_coeff = 9.0000e+06 2.000 6024.83
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_575 {
  chem_eq = "CH3 + C4H8_1 --> CH4 + C4H71_4"
  arrhenius_coeff = 9.0000e+04 2.000 6911.97
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_576 {
  chem_eq = "C2H3 + C4H8_1 --> C2H4 + C4H71_4"
  arrhenius_coeff = 2.0340e+05 2.000 5684.29
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_577 {
  chem_eq = "C2H5 + C4H8_1 --> C2H6 + C4H71_4"
  arrhenius_coeff = 6.0000e+04 2.000 9070.21
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_578 {
  chem_eq = "NC3H7 + C4H8_1 --> C3H8 + C4H71_4"
  arrhenius_coeff = 4.0500e+04 2.000 9070.21
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_579 {
  chem_eq = "IC3H7 + C4H8_1 --> C3H8 + C4H71_4"
  arrhenius_coeff = 4.0500e+04 2.000 12506.67
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_580 {
  chem_eq = "C3H5_S + C4H8_1 --> C3H6 + C4H71_4"
  arrhenius_coeff = 8.1000e+04 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_581 {
  chem_eq = "C3H5_T + C4H8_1 --> C3H6 + C4H71_4"
  arrhenius_coeff = 8.1000e+04 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_582 {
  chem_eq = "C3H5_A + C4H8_1 --> C3H6 + C4H71_4"
  arrhenius_coeff = 1.6170e+05 2.000 15810.16
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_583 {
  chem_eq = "C4H8_1 + C4H5 --> C4H71_4 + C4H6"
  arrhenius_coeff = 8.1000e+04 2.000 6342.26
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_584 {
  chem_eq = "C4H8_1 + C4H71_3 --> C4H8_2 + C4H71_4"
  arrhenius_coeff = 1.0200e+05 2.000 16763.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_585 {
  chem_eq = "IC4H7 + C4H8_1 --> IC4H8 + C4H71_4"
  arrhenius_coeff = 8.1000e+04 2.000 16763.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_586 {
  chem_eq = "C4H8_1 + C4H3 --> C4H71_4 + C4H4"
  arrhenius_coeff = 1.2270e+05 2.000 7129.99
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_587 {
  chem_eq = "C4H8_1 + C6H5 --> C4H71_4 + C6H6"
  arrhenius_coeff = 1.1100e+08 1.000 3309.53
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_588 {
  chem_eq = "C4H8_1 + C5H5 --> C4H71_4 + C5H6"
  arrhenius_coeff = 1.0200e+05 2.000 17080.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_589 {
  chem_eq = "C3H3 + C4H8_1 --> C3H4_A + C4H71_4"
  arrhenius_coeff = 8.1000e+04 2.000 16525.10
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_590 {
  chem_eq = "C4H8_1 + LC5H7 --> C4H71_4 + LC5H8"
  arrhenius_coeff = 8.1000e+04 2.000 18136.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_591 {
  chem_eq = "C4H8_1 + C7H7 --> C4H71_4 + C7H8"
  arrhenius_coeff = 4.0500e+04 2.000 19190.53
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_592 {
  chem_eq = "C4H8_1 + CH3C6H4 --> C4H71_4 + C7H8"
  arrhenius_coeff = 1.1100e+08 1.000 3864.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_593 {
  chem_eq = "H + C4H8_2 --> H2 + 0.03*CH3 + 0.03*C3H4_P + 0.97*C4H71_3"
  arrhenius_coeff = 1.8000e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_594 {
  chem_eq = "CH3 + C4H8_2 --> CH4 + 0.03*CH3 + 0.03*C3H4_P +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 1.8000e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_595 {
  chem_eq = "C2H3 + C4H8_2 --> 0.03*CH3 + C2H4 + 0.03*C3H4_P +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 4.0680e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_596 {
  chem_eq = "C2H5 + C4H8_2 --> 0.03*CH3 + C2H6 + 0.03*C3H4_P +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 1.2000e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_597 {
  chem_eq = "NC3H7 + C4H8_2 --> 0.03*CH3 + C3H8 + 0.03*C3H4_P +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 8.1000e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_598 {
  chem_eq = "IC3H7 + C4H8_2 --> 0.03*CH3 + C3H8 + 0.03*C3H4_P +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 8.1000e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_599 {
  chem_eq = "C3H5_S + C4H8_2 --> 0.03*CH3 + C3H6 + 0.03*C3H4_P +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 1.6200e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_600 {
  chem_eq = "C3H5_T + C4H8_2 --> 0.03*CH3 + C3H6 + 0.03*C3H4_P +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 1.6200e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_601 {
  chem_eq = "C3H5_A + C4H8_2 --> 0.03*CH3 + C3H6 + 0.03*C3H4_P +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 3.2340e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_602 {
  chem_eq = "C4H8_2 + C4H5 --> 0.03*CH3 + 0.03*C3H4_P + 0.97*C4H71_3 + C4H6"
  arrhenius_coeff = 1.6200e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_603 {
  chem_eq = "C4H8_2 + C4H71_4 --> 0.03*CH3 + 0.03*C3H4_P + C4H8_1 +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 8.1000e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_604 {
  chem_eq = "C4H8_2 + C4H71_3 --> 0.03*CH3 + 0.03*C3H4_P + C4H8_2 +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 2.0400e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_605 {
  chem_eq = "IC4H7 + C4H8_2 --> 0.03*CH3 + 0.03*C3H4_P + IC4H8 +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 1.6200e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_606 {
  chem_eq = "C4H8_2 + C4H3 --> 0.03*CH3 + 0.03*C3H4_P + 0.97*C4H71_3 + C4H4"
  arrhenius_coeff = 2.4540e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_607 {
  chem_eq = "C4H8_2 + C6H5 --> 0.03*CH3 + 0.03*C3H4_P + 0.97*C4H71_3 + C6H6"
  arrhenius_coeff = 2.2200e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_608 {
  chem_eq = "C4H8_2 + C5H5 --> 0.03*CH3 + 0.03*C3H4_P + 0.97*C4H71_3 + C5H6"
  arrhenius_coeff = 2.0400e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_609 {
  chem_eq = "C3H3 + C4H8_2 --> 0.03*CH3 + 0.03*C3H4_P + C3H4_A +" &
        " 0.97*C4H71_3"
  arrhenius_coeff = 1.6200e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_610 {
  chem_eq = "C4H8_2 + LC5H7 --> 0.03*CH3 + 0.03*C3H4_P + 0.97*C4H71_3 +" &
        " LC5H8"
  arrhenius_coeff = 1.6200e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_611 {
  chem_eq = "C4H8_2 + C7H7 --> 0.03*CH3 + 0.03*C3H4_P + 0.97*C4H71_3 + C7H8"
  arrhenius_coeff = 8.1000e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_612 {
  chem_eq = "C4H8_2 + CH3C6H4 --> 0.03*CH3 + 0.03*C3H4_P + 0.97*C4H71_3" &
        " + C7H8"
  arrhenius_coeff = 2.2200e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_613 {
  chem_eq = "H + IC4H8 --> H2 + IC4H7"
  arrhenius_coeff = 1.8000e+07 2.000 4389.88
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_614 {
  chem_eq = "CH3 + IC4H8 --> CH4 + IC4H7"
  arrhenius_coeff = 1.8000e+05 2.000 5638.84
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_615 {
  chem_eq = "C2H3 + IC4H8 --> C2H4 + IC4H7"
  arrhenius_coeff = 4.0680e+05 2.000 4567.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_616 {
  chem_eq = "C2H5 + IC4H8 --> C2H6 + IC4H7"
  arrhenius_coeff = 1.2000e+05 2.000 6963.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_617 {
  chem_eq = "NC3H7 + IC4H8 --> C3H8 + IC4H7"
  arrhenius_coeff = 8.1000e+04 2.000 6963.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_618 {
  chem_eq = "IC3H7 + IC4H8 --> C3H8 + IC4H7"
  arrhenius_coeff = 8.1000e+04 2.000 10163.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_619 {
  chem_eq = "IC4H8 + C4H5 --> IC4H7 + C4H6"
  arrhenius_coeff = 1.6200e+05 2.000 5167.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_620 {
  chem_eq = "IC4H8 + C4H71_4 --> IC4H7 + C4H8_1"
  arrhenius_coeff = 8.1000e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_621 {
  chem_eq = "IC4H8 + C4H71_3 --> IC4H7 + C4H8_2"
  arrhenius_coeff = 2.0400e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_622 {
  chem_eq = "IC4H8 + C4H3 --> IC4H7 + C4H4"
  arrhenius_coeff = 2.4540e+05 2.000 5838.84
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_623 {
  chem_eq = "IC4H8 + C6H5 --> IC4H7 + C6H6"
  arrhenius_coeff = 2.2200e+08 1.000 2667.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_624 {
  chem_eq = "IC4H8 + C5H5 --> IC4H7 + C5H6"
  arrhenius_coeff = 2.0400e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_625 {
  chem_eq = "C3H3 + IC4H8 --> C3H4_A + IC4H7"
  arrhenius_coeff = 1.6200e+05 2.000 14667.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_626 {
  chem_eq = "IC4H8 + LC5H7 --> IC4H7 + LC5H8"
  arrhenius_coeff = 1.6200e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_627 {
  chem_eq = "IC4H8 + C7H7 --> IC4H7 + C7H8"
  arrhenius_coeff = 8.1000e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_628 {
  chem_eq = "IC4H8 + CH3C6H4 --> IC4H7 + C7H8"
  arrhenius_coeff = 2.2200e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_629 {
  chem_eq = "C2H5 + C5H6 --> C2H6 + C5H5"
  arrhenius_coeff = 2.0000e+05 2.000 6980.29
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_630 {
  chem_eq = "NC3H7 + C5H6 --> C3H8 + C5H5"
  arrhenius_coeff = 1.3500e+05 2.000 6980.29
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_631 {
  chem_eq = "IC3H7 + C5H6 --> C3H8 + C5H5"
  arrhenius_coeff = 1.3500e+05 2.000 10180.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_632 {
  chem_eq = "C3H5_S + C5H6 --> C3H6 + C5H5"
  arrhenius_coeff = 2.7000e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_633 {
  chem_eq = "C3H5_T + C5H6 --> C3H6 + C5H5"
  arrhenius_coeff = 2.7000e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_634 {
  chem_eq = "C3H5_A + C5H6 --> C3H6 + C5H5"
  arrhenius_coeff = 5.3900e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_635 {
  chem_eq = "C4H5 + C5H6 --> C4H6 + C5H5"
  arrhenius_coeff = 2.7000e+05 2.000 5210.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_636 {
  chem_eq = "C4H71_4 + C5H6 --> C4H8_1 + C5H5"
  arrhenius_coeff = 1.3500e+05 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_637 {
  chem_eq = "C4H71_3 + C5H6 --> C4H8_2 + C5H5"
  arrhenius_coeff = 3.4000e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_638 {
  chem_eq = "IC4H7 + C5H6 --> IC4H8 + C5H5"
  arrhenius_coeff = 2.7000e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_639 {
  chem_eq = "C4H3 + C5H6 --> C4H4 + C5H5"
  arrhenius_coeff = 4.0900e+05 2.000 5880.69
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_640 {
  chem_eq = "C6H5 + C5H6 --> C6H6 + C5H5"
  arrhenius_coeff = 3.7000e+08 1.000 2710.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_641 {
  chem_eq = "C3H3 + C5H6 --> C3H4_A + C5H5"
  arrhenius_coeff = 2.7000e+05 2.000 14710.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_642 {
  chem_eq = "C5H6 + LC5H7 --> C5H5 + LC5H8"
  arrhenius_coeff = 2.7000e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_643 {
  chem_eq = "C5H6 + C7H7 --> C5H5 + C7H8"
  arrhenius_coeff = 1.3500e+05 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_644 {
  chem_eq = "C5H6 + CH3C6H4 --> C5H5 + C7H8"
  arrhenius_coeff = 3.7000e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_645 {
  chem_eq = "C2H3 + C5H5CH3 --> H + C2H4 + C6H6"
  arrhenius_coeff = 4.0680e+05 2.000 2730.26
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_646 {
  chem_eq = "C2H5 + C5H5CH3 --> H + C2H6 + C6H6"
  arrhenius_coeff = 1.2000e+05 2.000 5405.42
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_647 {
  chem_eq = "NC3H7 + C5H5CH3 --> H + C3H8 + C6H6"
  arrhenius_coeff = 8.1000e+04 2.000 5405.42
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_648 {
  chem_eq = "IC3H7 + C5H5CH3 --> H + C3H8 + C6H6"
  arrhenius_coeff = 8.1000e+04 2.000 8487.19
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_649 {
  chem_eq = "C3H5_S + C5H5CH3 --> H + C3H6 + C6H6"
  arrhenius_coeff = 1.6200e+05 2.000 3301.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_650 {
  chem_eq = "C3H5_T + C5H5CH3 --> H + C3H6 + C6H6"
  arrhenius_coeff = 1.6200e+05 2.000 3301.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_651 {
  chem_eq = "C3H5_A + C5H5CH3 --> H + C3H6 + C6H6"
  arrhenius_coeff = 3.2340e+05 2.000 11294.92
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_652 {
  chem_eq = "C4H5 + C5H5CH3 --> H + C4H6 + C6H6"
  arrhenius_coeff = 1.6200e+05 2.000 3301.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_653 {
  chem_eq = "C4H71_4 + C5H5CH3 --> H + C4H8_1 + C6H6"
  arrhenius_coeff = 8.1000e+04 2.000 5309.47
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_654 {
  chem_eq = "C4H71_3 + C5H5CH3 --> H + C4H8_2 + C6H6"
  arrhenius_coeff = 2.0400e+05 2.000 12168.47
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_655 {
  chem_eq = "IC4H7 + C5H5CH3 --> H + IC4H8 + C6H6"
  arrhenius_coeff = 1.6200e+05 2.000 12168.47
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_656 {
  chem_eq = "C4H3 + C5H5CH3 --> H + C4H4 + C6H6"
  arrhenius_coeff = 2.4540e+05 2.000 3969.14
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_657 {
  chem_eq = "C6H5 + C5H5CH3 --> H + 2.0*C6H6"
  arrhenius_coeff = 2.2200e+08 1.000 1067.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_658 {
  chem_eq = "C3H3 + C5H5CH3 --> H + C3H4_A + C6H6"
  arrhenius_coeff = 1.6200e+05 2.000 12459.86
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_659 {
  chem_eq = "C5H5CH3 + LC5H7 --> H + C6H6 + LC5H8"
  arrhenius_coeff = 1.6200e+05 2.000 13431.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_660 {
  chem_eq = "C5H5CH3 + C7H7 --> H + C6H6 + C7H8"
  arrhenius_coeff = 8.1000e+04 2.000 14404.74
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_661 {
  chem_eq = "C5H5CH3 + CH3C6H4 --> H + C6H6 + C7H8"
  arrhenius_coeff = 2.2200e+08 1.000 1067.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_662 {
  chem_eq = "C2H3 + C5H5CH3 --> H + C2H4 + FULVENE"
  arrhenius_coeff = 1.3560e+05 2.000 2730.26
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_663 {
  chem_eq = "C2H5 + C5H5CH3 --> H + C2H6 + FULVENE"
  arrhenius_coeff = 4.0000e+04 2.000 5405.42
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_664 {
  chem_eq = "NC3H7 + C5H5CH3 --> H + C3H8 + FULVENE"
  arrhenius_coeff = 2.7000e+04 2.000 5405.42
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_665 {
  chem_eq = "IC3H7 + C5H5CH3 --> H + C3H8 + FULVENE"
  arrhenius_coeff = 2.7000e+04 2.000 8487.19
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_666 {
  chem_eq = "C3H5_S + C5H5CH3 --> H + C3H6 + FULVENE"
  arrhenius_coeff = 5.4000e+04 2.000 3301.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_667 {
  chem_eq = "C3H5_T + C5H5CH3 --> H + C3H6 + FULVENE"
  arrhenius_coeff = 5.4000e+04 2.000 3301.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_668 {
  chem_eq = "C3H5_A + C5H5CH3 --> H + C3H6 + FULVENE"
  arrhenius_coeff = 1.0780e+05 2.000 11294.92
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_669 {
  chem_eq = "C4H5 + C5H5CH3 --> H + C4H6 + FULVENE"
  arrhenius_coeff = 5.4000e+04 2.000 3301.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_670 {
  chem_eq = "C4H71_4 + C5H5CH3 --> H + C4H8_1 + FULVENE"
  arrhenius_coeff = 2.7000e+04 2.000 5309.47
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_671 {
  chem_eq = "C4H71_3 + C5H5CH3 --> H + C4H8_2 + FULVENE"
  arrhenius_coeff = 6.8000e+04 2.000 12168.47
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_672 {
  chem_eq = "IC4H7 + C5H5CH3 --> H + IC4H8 + FULVENE"
  arrhenius_coeff = 5.4000e+04 2.000 12168.47
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_673 {
  chem_eq = "C4H3 + C5H5CH3 --> H + C4H4 + FULVENE"
  arrhenius_coeff = 8.1800e+04 2.000 3969.14
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_674 {
  chem_eq = "C6H5 + C5H5CH3 --> H + C6H6 + FULVENE"
  arrhenius_coeff = 7.4000e+07 1.000 1067.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_675 {
  chem_eq = "C3H3 + C5H5CH3 --> H + C3H4_A + FULVENE"
  arrhenius_coeff = 5.4000e+04 2.000 12459.86
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_676 {
  chem_eq = "C5H5CH3 + LC5H7 --> H + FULVENE + LC5H8"
  arrhenius_coeff = 5.4000e+04 2.000 13431.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_677 {
  chem_eq = "C5H5CH3 + C7H7 --> H + FULVENE + C7H8"
  arrhenius_coeff = 2.7000e+04 2.000 14404.74
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_678 {
  chem_eq = "C5H5CH3 + CH3C6H4 --> H + FULVENE + C7H8"
  arrhenius_coeff = 7.4000e+07 1.000 1067.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_679 {
  chem_eq = "C2H5 + C6H6 --> C2H6 + C6H5"
  arrhenius_coeff = 3.2000e+05 2.000 12330.09
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_680 {
  chem_eq = "NC3H7 + C6H6 --> C3H8 + C6H5"
  arrhenius_coeff = 2.1600e+05 2.000 12330.09
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_681 {
  chem_eq = "IC3H7 + C6H6 --> C3H8 + C6H5"
  arrhenius_coeff = 2.1600e+05 2.000 16129.68
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_682 {
  chem_eq = "C3H5_S + C6H6 --> C3H6 + C6H5"
  arrhenius_coeff = 4.3200e+05 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_683 {
  chem_eq = "C3H5_T + C6H6 --> C3H6 + C6H5"
  arrhenius_coeff = 4.3200e+05 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_684 {
  chem_eq = "C3H5_A + C6H6 --> C3H6 + C6H5"
  arrhenius_coeff = 8.6240e+05 2.000 20432.90
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_685 {
  chem_eq = "C4H5 + C6H6 --> C4H6 + C6H5"
  arrhenius_coeff = 4.3200e+05 2.000 8207.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_686 {
  chem_eq = "C4H71_4 + C6H6 --> C4H8_1 + C6H5"
  arrhenius_coeff = 2.1600e+05 2.000 13144.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_687 {
  chem_eq = "C4H71_3 + C6H6 --> C4H8_2 + C6H5"
  arrhenius_coeff = 5.4400e+05 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_688 {
  chem_eq = "IC4H7 + C6H6 --> IC4H8 + C6H5"
  arrhenius_coeff = 4.3200e+05 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_689 {
  chem_eq = "C4H3 + C6H6 --> C4H4 + C6H5"
  arrhenius_coeff = 6.5440e+05 2.000 9170.82
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_690 {
  chem_eq = "C6H6 + C5H5 --> C6H5 + C5H6"
  arrhenius_coeff = 5.4400e+05 2.000 21810.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_691 {
  chem_eq = "C3H3 + C6H6 --> C3H4_A + C6H5"
  arrhenius_coeff = 4.3200e+05 2.000 19439.11
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_692 {
  chem_eq = "C6H6 + LC5H7 --> C6H5 + LC5H8"
  arrhenius_coeff = 4.3200e+05 2.000 22952.95
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_693 {
  chem_eq = "C6H6 + C7H7 --> C6H5 + C7H8"
  arrhenius_coeff = 2.1600e+05 2.000 24090.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_694 {
  chem_eq = "C6H6 + CH3C6H4 --> C6H5 + C7H8"
  arrhenius_coeff = 5.9200e+08 1.000 6728.39
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_695 {
  chem_eq = "C3H6 + C3H3 --> CH3 + C5H6"
  arrhenius_coeff = 2.0000e+11 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_696 {
  chem_eq = "C3H3 + IC4H8 --> CH3 + C5H5CH3"
  arrhenius_coeff = 2.0000e+11 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_697 {
  chem_eq = "2.0*IC4H8 --> 2.0*C4H8_2"
  arrhenius_coeff = 3.0000e+16 0.000 65000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_697_rev {
  chem_eq = " 2.0*C4H8_2 --> 2.0*IC4H8 "
  arrhenius_coeff = 3.0000e+16 0.000 65000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_698 {
  chem_eq = "H + B1M2 --> CH3 + IC4H8"
  arrhenius_coeff = 5.0000e+13 0.000 1000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_698_rev {
  chem_eq = " CH3 + IC4H8 --> H + B1M2 "
  arrhenius_coeff = 5.0000e+13 0.000 1000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_699 {
  chem_eq = "IC4H8 + IC4H7 --> 2.0*H2 + CH3 + C7H8"
  arrhenius_coeff = 1.2000e+10 0.000 20500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_700 {
  chem_eq = "IC4H8 + IC4H7 --> C2H4 + IC3H7 + C3H4_A"
  arrhenius_coeff = 1.2000e+10 0.000 20500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_701 {
  chem_eq = "IC4H8 + IC4H7 --> H + C3H4_A + B1M2"
  arrhenius_coeff = 4.0000e+10 0.000 20500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_702 {
  chem_eq = "IC4H8 + IC4H7 --> 1.01*H2 + C2H5 + 0.01*C6H6 +" &
        " 0.99*C5H5CH3"
  arrhenius_coeff = 1.0000e+10 0.000 20500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_703 {
  chem_eq = "C3H6 + IC4H7 --> 1.2*H2 + 0.6*CH3 + 0.4*C2H5 + 0.6*C6H6 +" &
        " 0.4*LC5H8"
  arrhenius_coeff = 8.0000e+10 0.000 20500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_704 {
  chem_eq = "CH3 + C4H8_2 --> H + B1M3"
  arrhenius_coeff = 2.5000e+10 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_704_rev {
  chem_eq = " H + B1M3 --> CH3 + C4H8_2 "
  arrhenius_coeff = 2.5000e+10 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_705 {
  chem_eq = "H + B1M2 --> C2H5 + C3H6"
  arrhenius_coeff = 1.2500e+13 0.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_705_rev {
  chem_eq = " C2H5 + C3H6 --> H + B1M2 "
  arrhenius_coeff = 1.2500e+13 0.000 2000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_706 {
  chem_eq = "H + B1M3 --> C2H4 + IC3H7"
  arrhenius_coeff = 1.2500e+13 0.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_706_rev {
  chem_eq = " C2H4 + IC3H7 --> H + B1M3 "
  arrhenius_coeff = 1.2500e+13 0.000 2000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_707 {
  chem_eq = "C2H5 + C4H8_1 --> CH3 + B1M3"
  arrhenius_coeff = 7.5000e+10 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_707_rev {
  chem_eq = " CH3 + B1M3 --> C2H5 + C4H8_1 "
  arrhenius_coeff = 7.5000e+10 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_708 {
  chem_eq = "C2H5 + C4H8_1 --> 0.95*C2H4 + 0.05*NC3H7 + 0.05*C3H6 +" &
        " 0.8*PC4H9 + 0.15*SC4H9"
  arrhenius_coeff = 1.2000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_709 {
  chem_eq = "C2H5 + C4H8_2 --> C2H4 + SC4H9"
  arrhenius_coeff = 1.5000e+10 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_709_rev {
  chem_eq = " C2H4 + SC4H9 --> C2H5 + C4H8_2 "
  arrhenius_coeff = 1.5000e+10 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_710 {
  chem_eq = "C4H8_1 + C4H71_3 --> H2 + C2H5 + CYC6H8"
  arrhenius_coeff = 2.5000e+10 0.000 12500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_711 {
  chem_eq = "C4H8_2 + C4H71_3 --> H + C3H6 + LC5H8"
  arrhenius_coeff = 3.0000e+10 0.000 13500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_712 {
  chem_eq = "C4H8_2 + C4H71_3 --> H2 + C2H5 + C5H5CH3"
  arrhenius_coeff = 2.0000e+10 0.000 13500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_713 {
  chem_eq = "C4H8_2 + C4H71_3 --> H2 + C2H5 + CYC6H8"
  arrhenius_coeff = 2.5000e+10 0.000 13000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_714 {
  chem_eq = "2.0*C4H71_3 --> C4H8_1 + C4H6"
  arrhenius_coeff = 6.0000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_714_rev {
  chem_eq = " C4H8_1 + C4H6 --> 2.0*C4H71_3 "
  arrhenius_coeff = 6.0000e+11 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_715 {
  chem_eq = "C4H71_3 --> C2H4 + C2H3"
  arrhenius_coeff = 1.5000e+13 0.000 50000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_715_rev {
  chem_eq = " C2H4 + C2H3 --> C4H71_3 "
  arrhenius_coeff = 1.5000e+13 0.000 50000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_716 {
  chem_eq = "C3H6 + C4H71_3 --> 0.5*H2 + CH3 + 0.25*C2H4 + 0.5*C5H6 +" &
        " 0.5*CYC6H10"
  arrhenius_coeff = 1.5000e+11 0.000 13500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_717 {
  chem_eq = "CH3 + IC4H7 --> B1M2"
  arrhenius_coeff = 1.0000e+12 0.000 3000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_717_rev {
  chem_eq = " B1M2 --> CH3 + IC4H7 "
  arrhenius_coeff = 1.0000e+12 0.000 3000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_718 {
  chem_eq = "CH3 + C4H71_3 --> B1M3"
  arrhenius_coeff = 5.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_718_rev {
  chem_eq = " B1M3 --> CH3 + C4H71_3 "
  arrhenius_coeff = 5.0000e+12 0.000 0.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_719 {
  chem_eq = "IC3H7 + C3H6 --> CH3 + B1M3"
  arrhenius_coeff = 1.0000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_719_rev {
  chem_eq = " CH3 + B1M3 --> IC3H7 + C3H6 "
  arrhenius_coeff = 1.0000e+11 0.000 7600.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_720 {
  chem_eq = "CH3 + NC4H10 --> CH4 + PC4H9"
  arrhenius_coeff = 1.8000e+05 2.000 7443.70
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_721 {
  chem_eq = "C2H3 + NC4H10 --> C2H4 + PC4H9"
  arrhenius_coeff = 4.0680e+05 2.000 6239.48
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_722 {
  chem_eq = "C2H5 + NC4H10 --> C2H6 + PC4H9"
  arrhenius_coeff = 1.2000e+05 2.000 9289.16
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_723 {
  chem_eq = "NC3H7 + NC4H10 --> C3H8 + PC4H9"
  arrhenius_coeff = 8.1000e+04 2.000 9289.16
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_724 {
  chem_eq = "IC3H7 + NC4H10 --> C3H8 + PC4H9"
  arrhenius_coeff = 8.1000e+04 2.000 12725.62
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_725 {
  chem_eq = "C3H5_S + NC4H10 --> C3H6 + PC4H9"
  arrhenius_coeff = 1.6200e+05 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_726 {
  chem_eq = "C3H5_T + NC4H10 --> C3H6 + PC4H9"
  arrhenius_coeff = 1.6200e+05 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_727 {
  chem_eq = "C3H5_A + NC4H10 --> C3H6 + PC4H9"
  arrhenius_coeff = 3.2340e+05 2.000 15810.16
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_728 {
  chem_eq = "NC4H10 + C4H5 --> PC4H9 + C4H6"
  arrhenius_coeff = 1.6200e+05 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_729 {
  chem_eq = "NC4H10 + C4H71_4 --> PC4H9 + C4H8_1"
  arrhenius_coeff = 8.1000e+04 2.000 9181.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_730 {
  chem_eq = "NC4H10 + C4H71_3 --> PC4H9 + C4H8_2"
  arrhenius_coeff = 2.0400e+05 2.000 16763.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_731 {
  chem_eq = "NC4H10 + IC4H7 --> PC4H9 + IC4H8"
  arrhenius_coeff = 1.6200e+05 2.000 16763.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_732 {
  chem_eq = "NC4H10 + C4H3 --> PC4H9 + C4H4"
  arrhenius_coeff = 2.4540e+05 2.000 7661.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_733 {
  chem_eq = "NC4H10 + C6H5 --> PC4H9 + C6H6"
  arrhenius_coeff = 2.2200e+08 1.000 3864.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_734 {
  chem_eq = "NC4H10 + C5H5 --> PC4H9 + C5H6"
  arrhenius_coeff = 2.0400e+05 2.000 17080.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_735 {
  chem_eq = "C3H3 + NC4H10 --> C3H4_A + PC4H9"
  arrhenius_coeff = 1.6200e+05 2.000 17080.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_736 {
  chem_eq = "NC4H10 + LC5H7 --> PC4H9 + LC5H8"
  arrhenius_coeff = 1.6200e+05 2.000 18136.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_737 {
  chem_eq = "NC4H10 + C7H7 --> PC4H9 + C7H8"
  arrhenius_coeff = 8.1000e+04 2.000 19190.53
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_738 {
  chem_eq = "NC4H10 + CH3C6H4 --> PC4H9 + C7H8"
  arrhenius_coeff = 2.2200e+08 1.000 3864.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_739 {
  chem_eq = "CH3 + NC4H10 --> CH4 + SC4H9"
  arrhenius_coeff = 1.2000e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_740 {
  chem_eq = "C2H3 + NC4H10 --> C2H4 + SC4H9"
  arrhenius_coeff = 2.7120e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_741 {
  chem_eq = "C2H5 + NC4H10 --> C2H6 + SC4H9"
  arrhenius_coeff = 8.0000e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_742 {
  chem_eq = "NC3H7 + NC4H10 --> C3H8 + SC4H9"
  arrhenius_coeff = 5.4000e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_743 {
  chem_eq = "IC3H7 + NC4H10 --> C3H8 + SC4H9"
  arrhenius_coeff = 5.4000e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_744 {
  chem_eq = "C3H5_S + NC4H10 --> C3H6 + SC4H9"
  arrhenius_coeff = 1.0800e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_745 {
  chem_eq = "C3H5_T + NC4H10 --> C3H6 + SC4H9"
  arrhenius_coeff = 1.0800e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_746 {
  chem_eq = "C3H5_A + NC4H10 --> C3H6 + SC4H9"
  arrhenius_coeff = 2.1560e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_747 {
  chem_eq = "NC4H10 + C4H5 --> SC4H9 + C4H6"
  arrhenius_coeff = 1.0800e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_748 {
  chem_eq = "NC4H10 + C4H71_4 --> SC4H9 + C4H8_1"
  arrhenius_coeff = 5.4000e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_749 {
  chem_eq = "NC4H10 + C4H71_3 --> SC4H9 + C4H8_2"
  arrhenius_coeff = 1.3600e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_750 {
  chem_eq = "NC4H10 + IC4H7 --> SC4H9 + IC4H8"
  arrhenius_coeff = 1.0800e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_751 {
  chem_eq = "NC4H10 + C4H3 --> SC4H9 + C4H4"
  arrhenius_coeff = 1.6360e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_752 {
  chem_eq = "NC4H10 + C6H5 --> SC4H9 + C6H6"
  arrhenius_coeff = 1.4800e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_753 {
  chem_eq = "NC4H10 + C5H5 --> SC4H9 + C5H6"
  arrhenius_coeff = 1.3600e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_754 {
  chem_eq = "C3H3 + NC4H10 --> C3H4_A + SC4H9"
  arrhenius_coeff = 1.0800e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_755 {
  chem_eq = "NC4H10 + LC5H7 --> SC4H9 + LC5H8"
  arrhenius_coeff = 1.0800e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_756 {
  chem_eq = "NC4H10 + C7H7 --> SC4H9 + C7H8"
  arrhenius_coeff = 5.4000e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_757 {
  chem_eq = "NC4H10 + CH3C6H4 --> SC4H9 + C7H8"
  arrhenius_coeff = 1.4800e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_758 {
  chem_eq = "C2H3 + IC4H10 --> C2H4 + IC4H9"
  arrhenius_coeff = 6.1020e+05 2.000 6239.48
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_759 {
  chem_eq = "NC3H7 + IC4H10 --> C3H8 + IC4H9"
  arrhenius_coeff = 1.2150e+05 2.000 9289.16
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_760 {
  chem_eq = "IC3H7 + IC4H10 --> C3H8 + IC4H9"
  arrhenius_coeff = 1.2150e+05 2.000 12725.62
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_761 {
  chem_eq = "C3H5_S + IC4H10 --> C3H6 + IC4H9"
  arrhenius_coeff = 2.4300e+05 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_762 {
  chem_eq = "C3H5_T + IC4H10 --> C3H6 + IC4H9"
  arrhenius_coeff = 2.4300e+05 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_763 {
  chem_eq = "C3H5_A + IC4H10 --> C3H6 + IC4H9"
  arrhenius_coeff = 4.8510e+05 2.000 15810.16
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_764 {
  chem_eq = "IC4H10 + C4H5 --> IC4H9 + C4H6"
  arrhenius_coeff = 2.4300e+05 2.000 6897.45
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_765 {
  chem_eq = "IC4H10 + C4H71_4 --> IC4H9 + C4H8_1"
  arrhenius_coeff = 1.2150e+05 2.000 9181.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_766 {
  chem_eq = "IC4H10 + C4H71_3 --> IC4H9 + C4H8_2"
  arrhenius_coeff = 3.0600e+05 2.000 16763.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_767 {
  chem_eq = "IC4H10 + IC4H7 --> IC4H9 + IC4H8"
  arrhenius_coeff = 2.4300e+05 2.000 16763.05
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_768 {
  chem_eq = "IC4H10 + C4H3 --> IC4H9 + C4H4"
  arrhenius_coeff = 3.6810e+05 2.000 7661.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_769 {
  chem_eq = "IC4H10 + C6H5 --> IC4H9 + C6H6"
  arrhenius_coeff = 3.3300e+08 1.000 3864.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_770 {
  chem_eq = "IC4H10 + C5H5 --> IC4H9 + C5H6"
  arrhenius_coeff = 3.0600e+05 2.000 17080.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_771 {
  chem_eq = "C3H3 + IC4H10 --> C3H4_A + IC4H9"
  arrhenius_coeff = 2.4300e+05 2.000 17080.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_772 {
  chem_eq = "IC4H10 + LC5H7 --> IC4H9 + LC5H8"
  arrhenius_coeff = 2.4300e+05 2.000 18136.38
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_773 {
  chem_eq = "IC4H10 + C7H7 --> IC4H9 + C7H8"
  arrhenius_coeff = 1.2150e+05 2.000 19190.53
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_774 {
  chem_eq = "IC4H10 + CH3C6H4 --> IC4H9 + C7H8"
  arrhenius_coeff = 3.3300e+08 1.000 3864.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_775 {
  chem_eq = "H + B1M3 --> H2 + 0.5*H + 0.5*C2H3 + 0.5*C3H6 + 0.5*LC5H8"
  arrhenius_coeff = 3.0000e+06 2.000 2825.33
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_776 {
  chem_eq = "CH3 + B1M3 --> 0.5*H + CH4 + 0.5*C2H3 + 0.5*C3H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 3.0000e+04 2.000 3778.15
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_777 {
  chem_eq = "C2H3 + B1M3 --> 0.5*H + C2H4 + 0.5*C2H3 + 0.5*C3H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 6.7800e+04 2.000 2730.26
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_778 {
  chem_eq = "C2H5 + B1M3 --> 0.5*H + C2H6 + 0.5*C2H3 + 0.5*C3H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 2.0000e+04 2.000 5405.42
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_779 {
  chem_eq = "NC3H7 + B1M3 --> 0.5*H + 0.5*C2H3 + C3H8 + 0.5*C3H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.3500e+04 2.000 5405.42
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_780 {
  chem_eq = "IC3H7 + B1M3 --> 0.5*H + 0.5*C2H3 + C3H8 + 0.5*C3H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.3500e+04 2.000 8487.19
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_781 {
  chem_eq = "C3H5_S + B1M3 --> 0.5*H + 0.5*C2H3 + 1.5*C3H6 + 0.5*LC5H8"
  arrhenius_coeff = 2.7000e+04 2.000 3301.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_782 {
  chem_eq = "C3H5_T + B1M3 --> 0.5*H + 0.5*C2H3 + 1.5*C3H6 + 0.5*LC5H8"
  arrhenius_coeff = 2.7000e+04 2.000 3301.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_783 {
  chem_eq = "C3H5_A + B1M3 --> 0.5*H + 0.5*C2H3 + 1.5*C3H6 + 0.5*LC5H8"
  arrhenius_coeff = 5.3900e+04 2.000 11294.92
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_784 {
  chem_eq = "C4H5 + B1M3 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + C4H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 2.7000e+04 2.000 3301.28
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_785 {
  chem_eq = "C4H71_4 + B1M3 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + C4H8_1 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.3500e+04 2.000 5309.47
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_786 {
  chem_eq = "C4H71_3 + B1M3 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + C4H8_2 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 3.4000e+04 2.000 12168.47
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_787 {
  chem_eq = "IC4H7 + B1M3 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + IC4H8 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 2.7000e+04 2.000 12168.47
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_788 {
  chem_eq = "C4H3 + B1M3 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + C4H4 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 4.0900e+04 2.000 3969.14
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_789 {
  chem_eq = "C6H5 + B1M3 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + C6H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 3.7000e+07 1.000 1067.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_790 {
  chem_eq = "C5H5 + B1M3 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + C5H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 3.4000e+04 2.000 12459.86
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_791 {
  chem_eq = "C3H3 + B1M3 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + C3H4_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 2.7000e+04 2.000 12459.86
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_792 {
  chem_eq = "B1M3 + LC5H7 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + 1.5*LC5H8"
  arrhenius_coeff = 2.7000e+04 2.000 13431.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_793 {
  chem_eq = "B1M3 + C7H7 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + 0.5*LC5H8 + C7H8"
  arrhenius_coeff = 1.3500e+04 2.000 14404.74
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_794 {
  chem_eq = "B1M3 + CH3C6H4 --> 0.5*H + 0.5*C2H3 + 0.5*C3H6 + 0.5*LC5H8" &
        " + C7H8"
  arrhenius_coeff = 3.7000e+07 1.000 1067.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_795 {
  chem_eq = "H + B1M2 --> H2 + 0.5*H + 0.5*C2H4 + 0.5*C3H5_S +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 2.1000e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_796 {
  chem_eq = "CH3 + B1M2 --> 0.5*H + CH4 + 0.5*C2H4 + 0.5*C3H5_S +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 2.1000e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_797 {
  chem_eq = "C2H3 + B1M2 --> 0.5*H + 1.5*C2H4 + 0.5*C3H5_S + 0.5*LC5H8"
  arrhenius_coeff = 4.7460e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_798 {
  chem_eq = "C2H5 + B1M2 --> 0.5*H + C2H6 + 0.5*C2H4 + 0.5*C3H5_S +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.4000e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_799 {
  chem_eq = "NC3H7 + B1M2 --> 0.5*H + 0.5*C2H4 + C3H8 + 0.5*C3H5_S +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 9.4500e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_800 {
  chem_eq = "IC3H7 + B1M2 --> 0.5*H + 0.5*C2H4 + C3H8 + 0.5*C3H5_S +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 9.4500e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_801 {
  chem_eq = "C3H5_S + B1M2 --> 0.5*H + 0.5*C2H4 + C3H6 + 0.5*C3H5_S +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.8900e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_802 {
  chem_eq = "C3H5_T + B1M2 --> 0.5*H + 0.5*C2H4 + C3H6 + 0.5*C3H5_S +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.8900e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_803 {
  chem_eq = "C3H5_A + B1M2 --> 0.5*H + 0.5*C2H4 + C3H6 + 0.5*C3H5_S +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 3.7730e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_804 {
  chem_eq = "C4H5 + B1M2 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S + C4H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.8900e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_805 {
  chem_eq = "C4H71_4 + B1M2 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S + C4H8_1 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 9.4500e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_806 {
  chem_eq = "C4H71_3 + B1M2 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S + C4H8_2 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 2.3800e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_807 {
  chem_eq = "IC4H7 + B1M2 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S + IC4H8 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.8900e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_808 {
  chem_eq = "C4H3 + B1M2 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S + C4H4 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 2.8630e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_809 {
  chem_eq = "C6H5 + B1M2 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S + C6H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 2.5900e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_810 {
  chem_eq = "C5H5 + B1M2 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S + C5H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 2.3800e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_811 {
  chem_eq = "C3H3 + B1M2 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S + C3H4_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.8900e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_812 {
  chem_eq = "B1M2 + LC5H7 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S + 1.5*LC5H8"
  arrhenius_coeff = 1.8900e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_813 {
  chem_eq = "B1M2 + C7H7 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S + 0.5*LC5H8 + C7H8"
  arrhenius_coeff = 9.4500e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_814 {
  chem_eq = "B1M2 + CH3C6H4 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_S +" &
        " 0.5*LC5H8 + C7H8"
  arrhenius_coeff = 2.5900e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_815 {
  chem_eq = "H + B1M2 --> H2 + 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A +" &
        " 0.3*C4H6"
  arrhenius_coeff = 9.0000e+06 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_816 {
  chem_eq = "CH3 + B1M2 --> CH4 + 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A +" &
        " 0.3*C4H6"
  arrhenius_coeff = 9.0000e+04 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_817 {
  chem_eq = "C2H3 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + C2H4 + 0.7*C3H4_A +" &
        " 0.3*C4H6"
  arrhenius_coeff = 2.0340e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_818 {
  chem_eq = "C2H5 + B1M2 --> 0.3*CH3 + C2H6 + 0.7*C2H5 + 0.7*C3H4_A +" &
        " 0.3*C4H6"
  arrhenius_coeff = 6.0000e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_819 {
  chem_eq = "NC3H7 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + C3H8 + 0.7*C3H4_A +" &
        " 0.3*C4H6"
  arrhenius_coeff = 4.0500e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_820 {
  chem_eq = "IC3H7 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + C3H8 + 0.7*C3H4_A +" &
        " 0.3*C4H6"
  arrhenius_coeff = 4.0500e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_821 {
  chem_eq = "C3H5_S + B1M2 --> 0.3*CH3 + 0.7*C2H5 + C3H6 + 0.7*C3H4_A +" &
        " 0.3*C4H6"
  arrhenius_coeff = 8.1000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_822 {
  chem_eq = "C3H5_T + B1M2 --> 0.3*CH3 + 0.7*C2H5 + C3H6 + 0.7*C3H4_A +" &
        " 0.3*C4H6"
  arrhenius_coeff = 8.1000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_823 {
  chem_eq = "C3H5_A + B1M2 --> 0.3*CH3 + 0.7*C2H5 + C3H6 + 0.7*C3H4_A +" &
        " 0.3*C4H6"
  arrhenius_coeff = 1.6170e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_824 {
  chem_eq = "C4H5 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A + 1.3*C4H6"
  arrhenius_coeff = 8.1000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_825 {
  chem_eq = "C4H71_4 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A + C4H8_1" &
        " + 0.3*C4H6"
  arrhenius_coeff = 4.0500e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_826 {
  chem_eq = "C4H71_3 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A + C4H8_2" &
        " + 0.3*C4H6"
  arrhenius_coeff = 1.0200e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_827 {
  chem_eq = "IC4H7 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A + IC4H8 +" &
        " 0.3*C4H6"
  arrhenius_coeff = 8.1000e+04 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_828 {
  chem_eq = "C4H3 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A + 0.3*C4H6" &
        " + C4H4"
  arrhenius_coeff = 1.2270e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_829 {
  chem_eq = "C6H5 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A + 0.3*C4H6" &
        " + C6H6"
  arrhenius_coeff = 1.1100e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_830 {
  chem_eq = "C5H5 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A + 0.3*C4H6" &
        " + C5H6"
  arrhenius_coeff = 1.0200e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_831 {
  chem_eq = "C3H3 + B1M2 --> 0.3*CH3 + 0.7*C2H5 + 1.7*C3H4_A + 0.3*C4H6"
  arrhenius_coeff = 8.1000e+04 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_832 {
  chem_eq = "B1M2 + LC5H7 --> 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A + 0.3*C4H6" &
        " + LC5H8"
  arrhenius_coeff = 8.1000e+04 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_833 {
  chem_eq = "B1M2 + C7H7 --> 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A + 0.3*C4H6" &
        " + C7H8"
  arrhenius_coeff = 4.0500e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_834 {
  chem_eq = "B1M2 + CH3C6H4 --> 0.3*CH3 + 0.7*C2H5 + 0.7*C3H4_A +" &
        " 0.3*C4H6 + C7H8"
  arrhenius_coeff = 1.1100e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_835 {
  chem_eq = "NC5H11 --> 0.2*CH3 + 0.55*C2H5 + 0.25*C2H4 + 0.25*NC3H7 +" &
        " 0.55*C3H6 + 0.2*C4H8_1"
  arrhenius_coeff = 3.3000e+13 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_836 {
  chem_eq = "NC5H10 --> C2H5 + C3H5_A"
  arrhenius_coeff = 5.0000e+15 0.000 71000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_836_rev {
  chem_eq = " C2H5 + C3H5_A --> NC5H10 "
  arrhenius_coeff = 5.0000e+15 0.000 71000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_837 {
  chem_eq = "NC5H10 --> H + CH3 + C4H6"
  arrhenius_coeff = 6.0000e+18 0.000 93000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_838 {
  chem_eq = "NC5H10 --> CH3 + C4H71_3"
  arrhenius_coeff = 5.0000e+16 0.000 84000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_838_rev {
  chem_eq = " CH3 + C4H71_3 --> NC5H10 "
  arrhenius_coeff = 5.0000e+16 0.000 84000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_839 {
  chem_eq = "NC5H10 --> C2H3 + NC3H7"
  arrhenius_coeff = 1.0000e+17 0.000 99000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_839_rev {
  chem_eq = " C2H3 + NC3H7 --> NC5H10 "
  arrhenius_coeff = 1.0000e+17 0.000 99000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_840 {
  chem_eq = "NC5H10 --> C2H4 + C3H6"
  arrhenius_coeff = 6.5000e+11 0.000 53000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_840_rev {
  chem_eq = " C2H4 + C3H6 --> NC5H10 "
  arrhenius_coeff = 6.5000e+11 0.000 53000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_841 {
  chem_eq = "H + NC5H10 --> NC5H11"
  arrhenius_coeff = 2.0000e+13 0.000 2500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_841_rev {
  chem_eq = " NC5H11 --> H + NC5H10 "
  arrhenius_coeff = 2.0000e+13 0.000 2500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_842 {
  chem_eq = "H + NC5H10 --> H2 + CH3 + C4H6"
  arrhenius_coeff = 1.9800e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_843 {
  chem_eq = "CH3 + NC5H10 --> CH4 + CH3 + C4H6"
  arrhenius_coeff = 1.9800e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_844 {
  chem_eq = "C2H3 + NC5H10 --> CH3 + C2H4 + C4H6"
  arrhenius_coeff = 4.4748e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_845 {
  chem_eq = "C2H5 + NC5H10 --> CH3 + C2H6 + C4H6"
  arrhenius_coeff = 1.3200e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_846 {
  chem_eq = "NC3H7 + NC5H10 --> CH3 + C3H8 + C4H6"
  arrhenius_coeff = 8.9100e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_847 {
  chem_eq = "IC3H7 + NC5H10 --> CH3 + C3H8 + C4H6"
  arrhenius_coeff = 8.9100e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_848 {
  chem_eq = "C3H5_S + NC5H10 --> CH3 + C3H6 + C4H6"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_849 {
  chem_eq = "C3H5_T + NC5H10 --> CH3 + C3H6 + C4H6"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_850 {
  chem_eq = "C3H5_A + NC5H10 --> CH3 + C3H6 + C4H6"
  arrhenius_coeff = 3.5574e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_851 {
  chem_eq = "C4H5 + NC5H10 --> CH3 + 2.0*C4H6"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_852 {
  chem_eq = "C4H71_4 + NC5H10 --> CH3 + C4H8_1 + C4H6"
  arrhenius_coeff = 8.9100e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_853 {
  chem_eq = "C4H71_3 + NC5H10 --> CH3 + C4H8_2 + C4H6"
  arrhenius_coeff = 2.2440e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_854 {
  chem_eq = "IC4H7 + NC5H10 --> CH3 + IC4H8 + C4H6"
  arrhenius_coeff = 1.7820e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_855 {
  chem_eq = "C4H3 + NC5H10 --> CH3 + C4H6 + C4H4"
  arrhenius_coeff = 2.6994e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_856 {
  chem_eq = "C6H5 + NC5H10 --> CH3 + C4H6 + C6H6"
  arrhenius_coeff = 2.4420e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_857 {
  chem_eq = "C5H5 + NC5H10 --> CH3 + C4H6 + C5H6"
  arrhenius_coeff = 2.2440e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_858 {
  chem_eq = "C3H3 + NC5H10 --> CH3 + C3H4_A + C4H6"
  arrhenius_coeff = 1.7820e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_859 {
  chem_eq = "NC5H10 + LC5H7 --> CH3 + C4H6 + LC5H8"
  arrhenius_coeff = 1.7820e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_860 {
  chem_eq = "NC5H10 + C7H7 --> CH3 + C4H6 + C7H8"
  arrhenius_coeff = 8.9100e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_861 {
  chem_eq = "NC5H10 + CH3C6H4 --> CH3 + C4H6 + C7H8"
  arrhenius_coeff = 2.4420e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_862 {
  chem_eq = "H + NC5H10 --> H2 + 0.5*H + 0.5*C2H4 + 0.5*C3H5_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 7.2000e+06 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_863 {
  chem_eq = "CH3 + NC5H10 --> 0.5*H + CH4 + 0.5*C2H4 + 0.5*C3H5_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 7.2000e+04 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_864 {
  chem_eq = "C2H3 + NC5H10 --> 0.5*H + 1.5*C2H4 + 0.5*C3H5_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.6272e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_865 {
  chem_eq = "C2H5 + NC5H10 --> 0.5*H + C2H6 + 0.5*C2H4 + 0.5*C3H5_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 4.8000e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_866 {
  chem_eq = "NC3H7 + NC5H10 --> 0.5*H + 0.5*C2H4 + C3H8 + 0.5*C3H5_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 3.2400e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_867 {
  chem_eq = "IC3H7 + NC5H10 --> 0.5*H + 0.5*C2H4 + C3H8 + 0.5*C3H5_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 3.2400e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_868 {
  chem_eq = "C3H5_S + NC5H10 --> 0.5*H + 0.5*C2H4 + C3H6 + 0.5*C3H5_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 6.4800e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_869 {
  chem_eq = "C3H5_T + NC5H10 --> 0.5*H + 0.5*C2H4 + C3H6 + 0.5*C3H5_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 6.4800e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_870 {
  chem_eq = "C3H5_A + NC5H10 --> 0.5*H + 0.5*C2H4 + C3H6 + 0.5*C3H5_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 1.2936e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_871 {
  chem_eq = "C4H5 + NC5H10 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A + C4H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 6.4800e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_872 {
  chem_eq = "C4H71_4 + NC5H10 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A + C4H8_1" &
        " + 0.5*LC5H8"
  arrhenius_coeff = 3.2400e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_873 {
  chem_eq = "C4H71_3 + NC5H10 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A + C4H8_2" &
        " + 0.5*LC5H8"
  arrhenius_coeff = 8.1600e+04 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_874 {
  chem_eq = "IC4H7 + NC5H10 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A + IC4H8 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 6.4800e+04 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_875 {
  chem_eq = "C4H3 + NC5H10 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A + C4H4 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 9.8160e+04 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_876 {
  chem_eq = "C6H5 + NC5H10 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A + C6H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 8.8800e+07 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_877 {
  chem_eq = "C5H5 + NC5H10 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A + C5H6 +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 8.1600e+04 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_878 {
  chem_eq = "C3H3 + NC5H10 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A + C3H4_A +" &
        " 0.5*LC5H8"
  arrhenius_coeff = 6.4800e+04 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_879 {
  chem_eq = "NC5H10 + LC5H7 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A +" &
        " 1.5*LC5H8"
  arrhenius_coeff = 6.4800e+04 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_880 {
  chem_eq = "NC5H10 + C7H7 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A + 0.5*LC5H8" &
        " + C7H8"
  arrhenius_coeff = 3.2400e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_881 {
  chem_eq = "NC5H10 + CH3C6H4 --> 0.5*H + 0.5*C2H4 + 0.5*C3H5_A +" &
        " 0.5*LC5H8 + C7H8"
  arrhenius_coeff = 8.8800e+07 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_882 {
  chem_eq = "IC4H8 + C4H71_3 --> IC3H7 + LC5H8"
  arrhenius_coeff = 1.5000e+11 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_882_rev {
  chem_eq = " IC3H7 + LC5H8 --> IC4H8 + C4H71_3 "
  arrhenius_coeff = 1.5000e+11 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_883 {
  chem_eq = "C4H71_4 + C4H6 --> C3H5_A + LC5H8"
  arrhenius_coeff = 3.0000e+11 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_884 {
  chem_eq = "C3H6 + IC4H7 --> C2H5 + LC5H8"
  arrhenius_coeff = 2.5000e+11 0.000 18000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_884_rev {
  chem_eq = " C2H5 + LC5H8 --> C3H6 + IC4H7 "
  arrhenius_coeff = 2.5000e+11 0.000 18000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_885 {
  chem_eq = "C2H2 + C3H6 --> LC5H8"
  arrhenius_coeff = 1.5000e+11 0.000 32000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_885_rev {
  chem_eq = " LC5H8 --> C2H2 + C3H6 "
  arrhenius_coeff = 1.5000e+11 0.000 32000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_886 {
  chem_eq = "C2H4 + LC5H8 --> H + CH3 + CYC6H8"
  arrhenius_coeff = 2.5000e+10 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_887 {
  chem_eq = "C3H6 + C3H5_A --> CH3 + LC5H8"
  arrhenius_coeff = 2.0000e+10 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_887_rev {
  chem_eq = " CH3 + LC5H8 --> C3H6 + C3H5_A "
  arrhenius_coeff = 2.0000e+10 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_888 {
  chem_eq = "C3H5_A + C4H8_2 --> C2H5 + LC5H8"
  arrhenius_coeff = 7.0000e+10 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_888_rev {
  chem_eq = " C2H5 + LC5H8 --> C3H5_A + C4H8_2 "
  arrhenius_coeff = 7.0000e+10 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_889 {
  chem_eq = "LC5H8 --> H + LC5H7"
  arrhenius_coeff = 1.0000e+15 0.000 80000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_890 {
  chem_eq = "LC5H8 --> CH3 + C4H5"
  arrhenius_coeff = 1.0000e+16 0.000 87000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_891 {
  chem_eq = "LC5H8 --> C2H3 + C3H5_S"
  arrhenius_coeff = 7.5000e+16 0.000 105000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_892 {
  chem_eq = "C2H3 + C3H5_A --> 0.05*C2H4 + 0.05*C3H4_A + 0.95*LC5H8"
  arrhenius_coeff = 5.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_893 {
  chem_eq = "C6H6 + LC5H8 --> H2 + H + CH3 + C10H8"
  arrhenius_coeff = 5.0000e+11 0.000 44000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_894 {
  chem_eq = "H + LC5H8 --> 0.5*CH3 + 0.3*C2H4 + 0.2*C2H3 + 0.2*C3H6 +" &
        " 0.3*C3H5_S + 0.5*C4H6"
  arrhenius_coeff = 3.0000e+13 0.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_895 {
  chem_eq = "C2H3 + C3H6 --> H + LC5H8"
  arrhenius_coeff = 6.0000e+10 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_896 {
  chem_eq = "C2H3 + C4H8_1 --> CH3 + LC5H8"
  arrhenius_coeff = 5.0000e+11 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_897 {
  chem_eq = "C2H3 + C4H8_2 --> CH3 + LC5H8"
  arrhenius_coeff = 5.0000e+11 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_898 {
  chem_eq = "C2H3 + IC4H8 --> CH3 + LC5H8"
  arrhenius_coeff = 5.0000e+11 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_899 {
  chem_eq = "C2H3 + LC5H8 --> C2H5 + C5H6"
  arrhenius_coeff = 1.0000e+12 0.000 3000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_900 {
  chem_eq = "C2H4 + C3H5_A --> H + LC5H8"
  arrhenius_coeff = 3.3000e+10 0.000 22000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_900_rev {
  chem_eq = " H + LC5H8 --> C2H4 + C3H5_A "
  arrhenius_coeff = 3.3000e+10 0.000 22000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_901 {
  chem_eq = "C3H6 + C3H5_A --> 0.8*H2 + CH3 + 0.8*C5H6 + 0.2*LC5H8"
  arrhenius_coeff = 5.0000e+10 0.000 22000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_902 {
  chem_eq = "C3H5_A + IC4H8 --> C2H5 + LC5H8"
  arrhenius_coeff = 2.0000e+11 0.000 20000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_903 {
  chem_eq = "C3H6 + C3H5_S --> CH3 + LC5H8"
  arrhenius_coeff = 3.0000e+10 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_903_rev {
  chem_eq = " CH3 + LC5H8 --> C3H6 + C3H5_S "
  arrhenius_coeff = 3.0000e+10 0.000 6000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_904 {
  chem_eq = "H + LC5H8 --> H2 + 0.8*H + 0.2*C2H3 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6"
  arrhenius_coeff = 9.0000e+06 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_905 {
  chem_eq = "CH3 + LC5H8 --> 0.8*H + CH4 + 0.2*C2H3 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6"
  arrhenius_coeff = 9.0000e+04 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_906 {
  chem_eq = "C2H3 + LC5H8 --> 0.8*H + C2H4 + 0.2*C2H3 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6"
  arrhenius_coeff = 2.0340e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_907 {
  chem_eq = "C2H5 + LC5H8 --> 0.8*H + C2H6 + 0.2*C2H3 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6"
  arrhenius_coeff = 6.0000e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_908 {
  chem_eq = "NC3H7 + LC5H8 --> 0.8*H + 0.2*C2H3 + C3H8 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6"
  arrhenius_coeff = 4.0500e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_909 {
  chem_eq = "IC3H7 + LC5H8 --> 0.8*H + 0.2*C2H3 + C3H8 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6"
  arrhenius_coeff = 4.0500e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_910 {
  chem_eq = "C3H5_S + LC5H8 --> 0.8*H + 0.2*C2H3 + C3H6 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6"
  arrhenius_coeff = 8.1000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_911 {
  chem_eq = "C3H5_T + LC5H8 --> 0.8*H + 0.2*C2H3 + C3H6 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6"
  arrhenius_coeff = 8.1000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_912 {
  chem_eq = "C3H5_A + LC5H8 --> 0.8*H + 0.2*C2H3 + C3H6 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6"
  arrhenius_coeff = 1.6170e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_913 {
  chem_eq = "C4H5 + LC5H8 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P + 0.1*C3H4_A" &
        " + C4H6 + 0.8*C5H6"
  arrhenius_coeff = 8.1000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_914 {
  chem_eq = "C4H71_4 + LC5H8 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + C4H8_1 + 0.8*C5H6"
  arrhenius_coeff = 4.0500e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_915 {
  chem_eq = "C4H71_3 + LC5H8 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + C4H8_2 + 0.8*C5H6"
  arrhenius_coeff = 1.0200e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_916 {
  chem_eq = "IC4H7 + LC5H8 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + IC4H8 + 0.8*C5H6"
  arrhenius_coeff = 8.1000e+04 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_917 {
  chem_eq = "C4H3 + LC5H8 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P + 0.1*C3H4_A" &
        " + C4H4 + 0.8*C5H6"
  arrhenius_coeff = 1.2270e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_918 {
  chem_eq = "C6H5 + LC5H8 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P + 0.1*C3H4_A" &
        " + C6H6 + 0.8*C5H6"
  arrhenius_coeff = 1.1100e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_919 {
  chem_eq = "C5H5 + LC5H8 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P + 0.1*C3H4_A" &
        " + 1.8*C5H6"
  arrhenius_coeff = 1.0200e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_920 {
  chem_eq = "C3H3 + LC5H8 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P + 1.1*C3H4_A" &
        " + 0.8*C5H6"
  arrhenius_coeff = 8.1000e+04 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_921 {
  chem_eq = "LC5H8 + LC5H7 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6 + LC5H8"
  arrhenius_coeff = 8.1000e+04 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_922 {
  chem_eq = "LC5H8 + C7H7 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P + 0.1*C3H4_A" &
        " + 0.8*C5H6 + C7H8"
  arrhenius_coeff = 4.0500e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_923 {
  chem_eq = "LC5H8 + CH3C6H4 --> 0.8*H + 0.2*C2H3 + 0.1*C3H4_P +" &
        " 0.1*C3H4_A + 0.8*C5H6 + C7H8"
  arrhenius_coeff = 1.1100e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_924 {
  chem_eq = "H + C5H6 --> CYC5H7"
  arrhenius_coeff = 3.2800e+21 -3.140 1324.92
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 3.280000e+21" &
        " -3.140000e+00 1.325000e+03, 1.760000e+14 -6.800000e-01" &
        " -4.840000e+02, 6.650000e+10 6.400000e-01 -2.410000e+02"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_924_rev {
  chem_eq = " CYC5H7 --> H + C5H6 "
  arrhenius_coeff = 3.2800e+21 -3.140 1324.92
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 3.280000e+21" &
        " -3.140000e+00 1.325000e+03, 1.760000e+14 -6.800000e-01" &
        " -4.840000e+02, 6.650000e+10 6.400000e-01 -2.410000e+02"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_925 {
  chem_eq = "H + C5H6 --> CYC5H7"
  arrhenius_coeff = 1.9500e+69 -16.980 26774.85
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 1.950000e+69" &
        " -1.698000e+01 2.677500e+04, 6.500000e+72 -1.733000e+01" &
        " 3.553500e+04, 3.120000e+60 -1.322000e+01 3.626100e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_925_rev {
  chem_eq = " CYC5H7 --> H + C5H6 "
  arrhenius_coeff = 1.9500e+69 -16.980 26774.85
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 1.950000e+69" &
        " -1.698000e+01 2.677500e+04, 6.500000e+72 -1.733000e+01" &
        " 3.553500e+04, 3.120000e+60 -1.322000e+01 3.626100e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_926 {
  chem_eq = "H + C5H6 --> LC5H7"
  arrhenius_coeff = 8.0300e+42 -9.390 15954.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 8.030000e+42" &
        " -9.390000e+00 1.595400e+04, 2.580000e+21 -4.280000e+00" &
        " 2.290000e+03, 2.260000e+05 2.820000e+00 1.472300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_926_rev {
  chem_eq = " LC5H7 --> H + C5H6 "
  arrhenius_coeff = 8.0300e+42 -9.390 15954.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 8.030000e+42" &
        " -9.390000e+00 1.595400e+04, 2.580000e+21 -4.280000e+00" &
        " 2.290000e+03, 2.260000e+05 2.820000e+00 1.472300e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_927 {
  chem_eq = "H + C5H6 --> LC5H7"
  arrhenius_coeff = 5.3300e+68 -16.450 32251.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 arrhenius_press: 5.330000e+68 -1.645000e+01" &
        " 3.225100e+04, 5.000000e+99 -2.496000e+01 5.678800e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_927_rev {
  chem_eq = " LC5H7 --> H + C5H6 "
  arrhenius_coeff = 5.3300e+68 -16.450 32251.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 arrhenius_press: 5.330000e+68 -1.645000e+01" &
        " 3.225100e+04, 5.000000e+99 -2.496000e+01 5.678800e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_928 {
  chem_eq = "CYC5H7 --> LC5H7"
  arrhenius_coeff = 1.2100e+56 -13.260 60396.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 1.210000e+56" &
        " -1.326000e+01 6.039600e+04, 5.180000e+36 -7.010000e+00" &
        " 5.480000e+04, 8.070000e+22 -2.730000e+00 4.970500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_928_rev {
  chem_eq = " LC5H7 --> CYC5H7 "
  arrhenius_coeff = 1.2100e+56 -13.260 60396.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 1.210000e+56" &
        " -1.326000e+01 6.039600e+04, 5.180000e+36 -7.010000e+00" &
        " 5.480000e+04, 8.070000e+22 -2.730000e+00 4.970500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_929 {
  chem_eq = "CYC5H7 --> C2H4 + C3H3"
  arrhenius_coeff = 3.0100e+52 -12.760 88102.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 3.010000e+52" &
        " -1.276000e+01 8.810200e+04, 1.340000e+48 -1.010000e+01" &
        " 9.397800e+04, 1.060000e+27 -3.420000e+00 9.165400e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_930 {
  chem_eq = "C2H2 + C3H5_A --> CYC5H7"
  arrhenius_coeff = 8.5200e+53 -13.610 28381.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 8.520000e+53" &
        " -1.361000e+01 2.838100e+04, 1.040000e+34 -6.980000e+00" &
        " 2.182700e+04, 3.560000e+42 -9.160000e+00 2.795700e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_930_rev {
  chem_eq = " CYC5H7 --> C2H2 + C3H5_A "
  arrhenius_coeff = 8.5200e+53 -13.610 28381.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 8.520000e+53" &
        " -1.361000e+01 2.838100e+04, 1.040000e+34 -6.980000e+00" &
        " 2.182700e+04, 3.560000e+42 -9.160000e+00 2.795700e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_931 {
  chem_eq = "C2H2 + C3H5_A --> LC5H7"
  arrhenius_coeff = 1.4600e+53 -13.390 28100.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 1.460000e+53" &
        " -1.339000e+01 2.810000e+04, 2.850000e+41 -8.940000e+00" &
        " 2.811000e+04, 9.920000e+26 -4.410000e+00 2.322200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_931_rev {
  chem_eq = " LC5H7 --> C2H2 + C3H5_A "
  arrhenius_coeff = 1.4600e+53 -13.390 28100.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 1.460000e+53" &
        " -1.339000e+01 2.810000e+04, 2.850000e+41 -8.940000e+00" &
        " 2.811000e+04, 9.920000e+26 -4.410000e+00 2.322200e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_932 {
  chem_eq = "CH3 + C4H4 --> LC5H7"
  arrhenius_coeff = 7.1400e+98 -27.300 41665.48
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 7.140000e+98" &
        " -2.730000e+01 4.166548e+04, 9.530000e+77 -2.002000e+01" &
        " 3.909317e+04, 5.280000e+64 -1.545000e+01 3.770463e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_933 {
  chem_eq = "C2H3 + C3H4_A --> LC5H7"
  arrhenius_coeff = 4.0000e+52 -12.580 21620.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 4.000000e+52" &
        " -1.258000e+01 2.162000e+04, 1.130000e+32 -5.830000e+00" &
        " 1.639100e+04, 1.720000e+14 -4.000000e-01 9.033000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_933_rev {
  chem_eq = " LC5H7 --> C2H3 + C3H4_A "
  arrhenius_coeff = 4.0000e+52 -12.580 21620.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 4.000000e+52" &
        " -1.258000e+01 2.162000e+04, 1.130000e+32 -5.830000e+00" &
        " 1.639100e+04, 1.720000e+14 -4.000000e-01 9.033000e+03"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_934 {
  chem_eq = "C2H3 + C3H4_P --> LC5H7"
  arrhenius_coeff = 5.7500e+62 -16.800 21051.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-02" &
        " 1.000000e+00 1.000000e+02 arrhenius_press: 5.750000e+62" &
        " -1.680000e+01 2.105100e+04, 2.260000e+60 -1.500000e+01" &
        " 2.578800e+04, 3.150000e+57 -1.343000e+01 2.963500e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_935 {
  chem_eq = "H + CYC5H7 --> H2 + C5H6"
  arrhenius_coeff = 9.3000e+06 1.920 -743.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_935_rev {
  chem_eq = " H2 + C5H6 --> H + CYC5H7 "
  arrhenius_coeff = 9.3000e+06 1.920 -743.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_936 {
  chem_eq = "C5H5 + CYC5H7 --> 2.0*C5H6"
  arrhenius_coeff = 1.3000e+22 -2.720 3906.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_936_rev {
  chem_eq = " 2.0*C5H6 --> C5H5 + CYC5H7 "
  arrhenius_coeff = 1.3000e+22 -2.720 3906.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_937 {
  chem_eq = "C5H6 + LC5H7 --> C2H4 + C2H3 + C6H6"
  arrhenius_coeff = 1.0000e+12 0.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_938 {
  chem_eq = "DIALLYL --> 2.0*C3H5_A"
  arrhenius_coeff = 4.9000e+22 -2.060 63355.50
  press_rxn_param = PLOG "press_coeff: 1.000000e+00" &
        " 4.000000e+00 1.000000e+01 arrhenius_press: 5.070000e+47" &
        " -9.700000e+00 7.268000e+04, 4.220000e+39 -7.300000e+00" &
        " 6.939000e+04, 2.120000e+35 -6.000000e+00 6.762000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_938_rev {
  chem_eq = " 2.0*C3H5_A --> DIALLYL "
  arrhenius_coeff = 4.9000e+22 -2.060 63355.50
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e+00" &
        " 4.000000e+00 1.000000e+01 arrhenius_press: 5.070000e+47" &
        " -9.700000e+00 7.268000e+04, 4.220000e+39 -7.300000e+00" &
        " 6.939000e+04, 2.120000e+35 -6.000000e+00 6.762000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_939 {
  chem_eq = "C3H6 + C3H5_A --> H + DIALLYL"
  arrhenius_coeff = 2.0000e+10 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_940 {
  chem_eq = "RC6H9A --> RCYC6H9"
  arrhenius_coeff = 1.0000e+10 0.000 18000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_940_rev {
  chem_eq = " RCYC6H9 --> RC6H9A "
  arrhenius_coeff = 1.0000e+10 0.000 18000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_941 {
  chem_eq = "C3H5_A + C3H4_P --> H + CYC6H8"
  arrhenius_coeff = 1.5000e+11 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_941_rev {
  chem_eq = " H + CYC6H8 --> C3H5_A + C3H4_P "
  arrhenius_coeff = 1.5000e+11 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_942 {
  chem_eq = "C3H5_A + C3H4_A --> H + CYC6H8"
  arrhenius_coeff = 1.5000e+11 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_942_rev {
  chem_eq = " H + CYC6H8 --> C3H5_A + C3H4_A "
  arrhenius_coeff = 1.5000e+11 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_943 {
  chem_eq = "C3H6 + C3H5_A --> H + CYC6H10"
  arrhenius_coeff = 2.5000e+10 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_943_rev {
  chem_eq = " H + CYC6H10 --> C3H6 + C3H5_A "
  arrhenius_coeff = 2.5000e+10 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_944 {
  chem_eq = "C3H5_A + C4H8_1 --> CH3 + CYC6H10"
  arrhenius_coeff = 1.0000e+10 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_944_rev {
  chem_eq = " CH3 + CYC6H10 --> C3H5_A + C4H8_1 "
  arrhenius_coeff = 1.0000e+10 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_945 {
  chem_eq = "C3H5_A + C4H6 --> CH3 + CYC6H8"
  arrhenius_coeff = 1.0000e+11 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_945_rev {
  chem_eq = " CH3 + CYC6H8 --> C3H5_A + C4H6 "
  arrhenius_coeff = 1.0000e+11 0.000 15000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_946 {
  chem_eq = "C3H5_A + C3H5_S --> CYC6H10"
  arrhenius_coeff = 6.0000e+12 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_947 {
  chem_eq = "C3H6 + IC4H7 --> CH3 + CYC6H10"
  arrhenius_coeff = 3.5000e+11 0.000 18000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_947_rev {
  chem_eq = " CH3 + CYC6H10 --> C3H6 + IC4H7 "
  arrhenius_coeff = 3.5000e+11 0.000 18000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_948 {
  chem_eq = "H + DIALLYL --> H2 + RC6H9A"
  arrhenius_coeff = 3.6000e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_949 {
  chem_eq = "CH3 + DIALLYL --> CH4 + RC6H9A"
  arrhenius_coeff = 3.6000e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_950 {
  chem_eq = "C2H3 + DIALLYL --> C2H4 + RC6H9A"
  arrhenius_coeff = 8.1360e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_951 {
  chem_eq = "C2H5 + DIALLYL --> C2H6 + RC6H9A"
  arrhenius_coeff = 2.4000e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_952 {
  chem_eq = "NC3H7 + DIALLYL --> C3H8 + RC6H9A"
  arrhenius_coeff = 1.6200e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_953 {
  chem_eq = "IC3H7 + DIALLYL --> C3H8 + RC6H9A"
  arrhenius_coeff = 1.6200e+05 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_954 {
  chem_eq = "C3H5_S + DIALLYL --> C3H6 + RC6H9A"
  arrhenius_coeff = 3.2400e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_955 {
  chem_eq = "C3H5_T + DIALLYL --> C3H6 + RC6H9A"
  arrhenius_coeff = 3.2400e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_956 {
  chem_eq = "C3H5_A + DIALLYL --> C3H6 + RC6H9A"
  arrhenius_coeff = 6.4680e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_957 {
  chem_eq = "C4H5 + DIALLYL --> C4H6 + RC6H9A"
  arrhenius_coeff = 3.2400e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_958 {
  chem_eq = "C4H71_4 + DIALLYL --> C4H8_1 + RC6H9A"
  arrhenius_coeff = 1.6200e+05 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_959 {
  chem_eq = "C4H71_3 + DIALLYL --> C4H8_2 + RC6H9A"
  arrhenius_coeff = 4.0800e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_960 {
  chem_eq = "IC4H7 + DIALLYL --> IC4H8 + RC6H9A"
  arrhenius_coeff = 3.2400e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_961 {
  chem_eq = "C4H3 + DIALLYL --> C4H4 + RC6H9A"
  arrhenius_coeff = 4.9080e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_962 {
  chem_eq = "C6H5 + DIALLYL --> C6H6 + RC6H9A"
  arrhenius_coeff = 4.4400e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_963 {
  chem_eq = "C5H5 + DIALLYL --> C5H6 + RC6H9A"
  arrhenius_coeff = 4.0800e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_964 {
  chem_eq = "C3H3 + DIALLYL --> C3H4_A + RC6H9A"
  arrhenius_coeff = 3.2400e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_965 {
  chem_eq = "LC5H7 + DIALLYL --> LC5H8 + RC6H9A"
  arrhenius_coeff = 3.2400e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_966 {
  chem_eq = "DIALLYL + C7H7 --> RC6H9A + C7H8"
  arrhenius_coeff = 1.6200e+05 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_967 {
  chem_eq = "DIALLYL + CH3C6H4 --> RC6H9A + C7H8"
  arrhenius_coeff = 4.4400e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_968 {
  chem_eq = "CYC6H10 --> CH3 + LC5H7"
  arrhenius_coeff = 2.5000e+16 0.000 71687.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_968_rev {
  chem_eq = " CH3 + LC5H7 --> CYC6H10 "
  arrhenius_coeff = 2.5000e+16 0.000 71687.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_969 {
  chem_eq = "CYC6H10 --> C2H4 + C4H6"
  arrhenius_coeff = 1.0000e+15 0.000 67000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_969_rev {
  chem_eq = " C2H4 + C4H6 --> CYC6H10 "
  arrhenius_coeff = 1.0000e+15 0.000 67000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_970 {
  chem_eq = "CYC6H8 --> H2 + C6H6"
  arrhenius_coeff = 1.0000e+14 0.000 70000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_970_rev {
  chem_eq = " H2 + C6H6 --> CYC6H8 "
  arrhenius_coeff = 1.0000e+14 0.000 70000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_971 {
  chem_eq = "CYC6H10 --> H2 + CYC6H8"
  arrhenius_coeff = 1.0000e+14 0.000 69000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_971_rev {
  chem_eq = " H2 + CYC6H8 --> CYC6H10 "
  arrhenius_coeff = 1.0000e+14 0.000 69000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_972 {
  chem_eq = "C2H4 + C4H4 --> CYC6H8"
  arrhenius_coeff = 3.0000e+11 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_973 {
  chem_eq = "RCYC6H9 --> H + CYC6H8"
  arrhenius_coeff = 3.0000e+13 0.000 36000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_974 {
  chem_eq = "NC7H14 --> C3H5_A + PC4H9"
  arrhenius_coeff = 1.5000e+16 0.000 71000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_974_rev {
  chem_eq = " C3H5_A + PC4H9 --> NC7H14 "
  arrhenius_coeff = 1.5000e+16 0.000 71000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_975 {
  chem_eq = "NC7H14 --> NC3H7 + C4H71_4"
  arrhenius_coeff = 5.0000e+16 0.000 84000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_975_rev {
  chem_eq = " NC3H7 + C4H71_4 --> NC7H14 "
  arrhenius_coeff = 5.0000e+16 0.000 84000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_976 {
  chem_eq = "NC7H14 --> C2H5 + C2H4 + C3H5_A"
  arrhenius_coeff = 5.0000e+16 0.000 84500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_977 {
  chem_eq = "NC7H14 --> C3H6 + C4H8_1"
  arrhenius_coeff = 6.5000e+11 0.000 53000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_977_rev {
  chem_eq = " C3H6 + C4H8_1 --> NC7H14 "
  arrhenius_coeff = 6.5000e+11 0.000 53000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_978 {
  chem_eq = "NC7H14 --> C2H4 + NC5H10"
  arrhenius_coeff = 2.0000e+11 0.000 53000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_978_rev {
  chem_eq = " C2H4 + NC5H10 --> NC7H14 "
  arrhenius_coeff = 2.0000e+11 0.000 53000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_979 {
  chem_eq = "H + NC7H14 --> NC7H15"
  arrhenius_coeff = 2.5000e+13 0.000 2500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_979_rev {
  chem_eq = " NC7H15 --> H + NC7H14 "
  arrhenius_coeff = 2.5000e+13 0.000 2500.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_980 {
  chem_eq = "NC7H15 --> C2H4 + NC5H11"
  arrhenius_coeff = 6.2000e+12 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_980_rev {
  chem_eq = " C2H4 + NC5H11 --> NC7H15 "
  arrhenius_coeff = 6.2000e+12 0.000 30000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_981 {
  chem_eq = "NC7H15 --> C3H6 + PC4H9"
  arrhenius_coeff = 1.4200e+13 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_981_rev {
  chem_eq = " C3H6 + PC4H9 --> NC7H15 "
  arrhenius_coeff = 1.4200e+13 0.000 30000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_982 {
  chem_eq = "NC7H15 --> NC3H7 + C4H8_1"
  arrhenius_coeff = 7.5000e+12 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_982_rev {
  chem_eq = " NC3H7 + C4H8_1 --> NC7H15 "
  arrhenius_coeff = 7.5000e+12 0.000 30000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_983 {
  chem_eq = "NC7H15 --> C2H5 + NC5H10"
  arrhenius_coeff = 5.6000e+12 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_983_rev {
  chem_eq = " C2H5 + NC5H10 --> NC7H15 "
  arrhenius_coeff = 5.6000e+12 0.000 30000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_984 {
  chem_eq = "DIMEPTD --> H + CH3 + C5H5CH3"
  arrhenius_coeff = 1.0000e+16 0.000 76000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_985 {
  chem_eq = "C4H8_1 + C4H71_3 --> CH3 + DIMEPTD"
  arrhenius_coeff = 5.0000e+10 0.000 12500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_986 {
  chem_eq = "H + DIMEPTD --> 0.9*CH3 + 0.1*C3H5_S + 0.325*IC4H8 +" &
        " 0.9*LC5H8"
  arrhenius_coeff = 2.5000e+13 0.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_987 {
  chem_eq = "IC4H8 + IC4H7 --> CH3 + DIMEPTD"
  arrhenius_coeff = 2.5000e+10 0.000 20000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_988 {
  chem_eq = "H + NC7H14 --> H2 + NC3H7 + C4H6"
  arrhenius_coeff = 1.9800e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_989 {
  chem_eq = "CH3 + NC7H14 --> CH4 + NC3H7 + C4H6"
  arrhenius_coeff = 1.9800e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_990 {
  chem_eq = "C2H3 + NC7H14 --> C2H4 + NC3H7 + C4H6"
  arrhenius_coeff = 4.4748e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_991 {
  chem_eq = "C2H5 + NC7H14 --> C2H6 + NC3H7 + C4H6"
  arrhenius_coeff = 1.3200e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_992 {
  chem_eq = "NC3H7 + NC7H14 --> C3H8 + NC3H7 + C4H6"
  arrhenius_coeff = 8.9100e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_993 {
  chem_eq = "IC3H7 + NC7H14 --> C3H8 + NC3H7 + C4H6"
  arrhenius_coeff = 8.9100e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_994 {
  chem_eq = "C3H5_S + NC7H14 --> NC3H7 + C3H6 + C4H6"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_995 {
  chem_eq = "C3H5_T + NC7H14 --> NC3H7 + C3H6 + C4H6"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_996 {
  chem_eq = "C3H5_A + NC7H14 --> NC3H7 + C3H6 + C4H6"
  arrhenius_coeff = 3.5574e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_997 {
  chem_eq = "C4H5 + NC7H14 --> NC3H7 + 2.0*C4H6"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_998 {
  chem_eq = "C4H71_4 + NC7H14 --> NC3H7 + C4H8_1 + C4H6"
  arrhenius_coeff = 8.9100e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_999 {
  chem_eq = "C4H71_3 + NC7H14 --> NC3H7 + C4H8_2 + C4H6"
  arrhenius_coeff = 2.2440e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1000 {
  chem_eq = "IC4H7 + NC7H14 --> NC3H7 + IC4H8 + C4H6"
  arrhenius_coeff = 1.7820e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1001 {
  chem_eq = "C4H3 + NC7H14 --> NC3H7 + C4H6 + C4H4"
  arrhenius_coeff = 2.6994e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1002 {
  chem_eq = "C6H5 + NC7H14 --> NC3H7 + C4H6 + C6H6"
  arrhenius_coeff = 2.4420e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1003 {
  chem_eq = "C5H5 + NC7H14 --> NC3H7 + C4H6 + C5H6"
  arrhenius_coeff = 2.2440e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1004 {
  chem_eq = "C3H3 + NC7H14 --> NC3H7 + C3H4_A + C4H6"
  arrhenius_coeff = 1.7820e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1005 {
  chem_eq = "LC5H7 + NC7H14 --> NC3H7 + C4H6 + LC5H8"
  arrhenius_coeff = 1.7820e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1006 {
  chem_eq = "NC7H14 + C7H7 --> NC3H7 + C4H6 + C7H8"
  arrhenius_coeff = 8.9100e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1007 {
  chem_eq = "NC7H14 + CH3C6H4 --> NC3H7 + C4H6 + C7H8"
  arrhenius_coeff = 2.4420e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1008 {
  chem_eq = "H + NC7H14 --> H2 + 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 1.9200e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1009 {
  chem_eq = "CH3 + NC7H14 --> CH4 + 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 1.9200e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1010 {
  chem_eq = "C2H3 + NC7H14 --> 0.5*C2H5 + C2H4 + 0.25*C3H6 + 0.25*C3H5_A" &
        " + 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 4.3392e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1011 {
  chem_eq = "C2H5 + NC7H14 --> C2H6 + 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A" &
        " + 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 1.2800e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1012 {
  chem_eq = "NC3H7 + NC7H14 --> 0.5*C2H5 + C3H8 + 0.25*C3H6 +" &
        " 0.25*C3H5_A + 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 8.6400e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1013 {
  chem_eq = "IC3H7 + NC7H14 --> 0.5*C2H5 + C3H8 + 0.25*C3H6 +" &
        " 0.25*C3H5_A + 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 8.6400e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1014 {
  chem_eq = "C3H5_S + NC7H14 --> 0.5*C2H5 + 1.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 1.7280e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1015 {
  chem_eq = "C3H5_T + NC7H14 --> 0.5*C2H5 + 1.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 1.7280e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1016 {
  chem_eq = "C3H5_A + NC7H14 --> 0.5*C2H5 + 1.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 3.4496e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1017 {
  chem_eq = "C4H5 + NC7H14 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + C4H6 + 0.5*LC5H8"
  arrhenius_coeff = 1.7280e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1018 {
  chem_eq = "C4H71_4 + NC7H14 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 1.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 8.6400e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1019 {
  chem_eq = "C4H71_3 + NC7H14 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + C4H8_2 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 2.1760e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1020 {
  chem_eq = "IC4H7 + NC7H14 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " IC4H8 + 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 1.7280e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1021 {
  chem_eq = "C4H3 + NC7H14 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + C4H4 + 0.5*LC5H8"
  arrhenius_coeff = 2.6176e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1022 {
  chem_eq = "C6H5 + NC7H14 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + C6H6 + 0.5*LC5H8"
  arrhenius_coeff = 2.3680e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1023 {
  chem_eq = "C5H5 + NC7H14 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + C5H6 + 0.5*LC5H8"
  arrhenius_coeff = 2.1760e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1024 {
  chem_eq = "C3H3 + NC7H14 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " C3H4_A + 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8"
  arrhenius_coeff = 1.7280e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1025 {
  chem_eq = "LC5H7 + NC7H14 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + 1.5*LC5H8"
  arrhenius_coeff = 1.7280e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1026 {
  chem_eq = "NC7H14 + C7H7 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8 + C7H8"
  arrhenius_coeff = 8.6400e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1027 {
  chem_eq = "NC7H14 + CH3C6H4 --> 0.5*C2H5 + 0.25*C3H6 + 0.25*C3H5_A +" &
        " 0.25*C4H8_1 + 0.25*C4H71_4 + 0.5*LC5H8 + C7H8"
  arrhenius_coeff = 2.3680e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1028 {
  chem_eq = "H + DIMEPTD --> 1.8*H2 + 0.8*CH3 + 0.2*C2H3 + 0.8*C6H6 +" &
        " 0.2*LC5H8"
  arrhenius_coeff = 1.6200e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1029 {
  chem_eq = "CH3 + DIMEPTD --> 0.8*H2 + CH4 + 0.8*CH3 + 0.2*C2H3 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 1.6200e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1030 {
  chem_eq = "C2H3 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + C2H4 + 0.2*C2H3 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 3.6612e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1031 {
  chem_eq = "C2H5 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + C2H6 + 0.2*C2H3 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 1.0800e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1032 {
  chem_eq = "NC3H7 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + C3H8 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 7.2900e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1033 {
  chem_eq = "IC3H7 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + C3H8 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 7.2900e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1034 {
  chem_eq = "C3H5_S + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + C3H6 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 1.4580e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1035 {
  chem_eq = "C3H5_T + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + C3H6 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 1.4580e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1036 {
  chem_eq = "C3H5_A + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + C3H6 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 2.9106e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1037 {
  chem_eq = "C4H5 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + C4H6 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 1.4580e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1038 {
  chem_eq = "C4H71_4 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + C4H8_1" &
        " + 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 7.2900e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1039 {
  chem_eq = "C4H71_3 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + C4H8_2" &
        " + 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 1.8360e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1040 {
  chem_eq = "IC4H7 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + IC4H8 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 1.4580e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1041 {
  chem_eq = "C4H3 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + C4H4 +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 2.2086e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1042 {
  chem_eq = "C6H5 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + 1.8*C6H6 +" &
        " 0.2*LC5H8"
  arrhenius_coeff = 1.9980e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1043 {
  chem_eq = "C5H5 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + 0.8*C6H6 +" &
        " C5H6 + 0.2*LC5H8"
  arrhenius_coeff = 1.8360e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1044 {
  chem_eq = "C3H3 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + C3H4_A +" &
        " 0.8*C6H6 + 0.2*LC5H8"
  arrhenius_coeff = 1.4580e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1045 {
  chem_eq = "LC5H7 + DIMEPTD --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + 0.8*C6H6" &
        " + 1.2*LC5H8"
  arrhenius_coeff = 1.4580e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1046 {
  chem_eq = "DIMEPTD + C7H7 --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 + 0.8*C6H6 +" &
        " 0.2*LC5H8 + C7H8"
  arrhenius_coeff = 7.2900e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1047 {
  chem_eq = "DIMEPTD + CH3C6H4 --> 0.8*H2 + 0.8*CH3 + 0.2*C2H3 +" &
        " 0.8*C6H6 + 0.2*LC5H8 + C7H8"
  arrhenius_coeff = 1.9980e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1048 {
  chem_eq = "H + DIMEPTD --> H2 + 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7 +" &
        " 0.5*C4H6"
  arrhenius_coeff = 1.0800e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1049 {
  chem_eq = "CH3 + DIMEPTD --> CH4 + 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7" &
        " + 0.5*C4H6"
  arrhenius_coeff = 1.0800e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1050 {
  chem_eq = "C2H3 + DIMEPTD --> C2H4 + 0.5*C3H5_A + 0.5*C3H4_A +" &
        " 0.5*IC4H7 + 0.5*C4H6"
  arrhenius_coeff = 2.4408e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1051 {
  chem_eq = "C2H5 + DIMEPTD --> C2H6 + 0.5*C3H5_A + 0.5*C3H4_A +" &
        " 0.5*IC4H7 + 0.5*C4H6"
  arrhenius_coeff = 7.2000e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1052 {
  chem_eq = "NC3H7 + DIMEPTD --> C3H8 + 0.5*C3H5_A + 0.5*C3H4_A +" &
        " 0.5*IC4H7 + 0.5*C4H6"
  arrhenius_coeff = 4.8600e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1053 {
  chem_eq = "IC3H7 + DIMEPTD --> C3H8 + 0.5*C3H5_A + 0.5*C3H4_A +" &
        " 0.5*IC4H7 + 0.5*C4H6"
  arrhenius_coeff = 4.8600e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1054 {
  chem_eq = "C3H5_S + DIMEPTD --> C3H6 + 0.5*C3H5_A + 0.5*C3H4_A +" &
        " 0.5*IC4H7 + 0.5*C4H6"
  arrhenius_coeff = 9.7200e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1055 {
  chem_eq = "C3H5_T + DIMEPTD --> C3H6 + 0.5*C3H5_A + 0.5*C3H4_A +" &
        " 0.5*IC4H7 + 0.5*C4H6"
  arrhenius_coeff = 9.7200e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1056 {
  chem_eq = "C3H5_A + DIMEPTD --> C3H6 + 0.5*C3H5_A + 0.5*C3H4_A +" &
        " 0.5*IC4H7 + 0.5*C4H6"
  arrhenius_coeff = 1.9404e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1057 {
  chem_eq = "C4H5 + DIMEPTD --> 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7 +" &
        " 1.5*C4H6"
  arrhenius_coeff = 9.7200e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1058 {
  chem_eq = "C4H71_4 + DIMEPTD --> 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7 +" &
        " C4H8_1 + 0.5*C4H6"
  arrhenius_coeff = 4.8600e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1059 {
  chem_eq = "C4H71_3 + DIMEPTD --> 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7 +" &
        " C4H8_2 + 0.5*C4H6"
  arrhenius_coeff = 1.2240e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1060 {
  chem_eq = "IC4H7 + DIMEPTD --> 0.5*C3H5_A + 0.5*C3H4_A + IC4H8 +" &
        " 0.5*IC4H7 + 0.5*C4H6"
  arrhenius_coeff = 9.7200e+04 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1061 {
  chem_eq = "C4H3 + DIMEPTD --> 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7 +" &
        " 0.5*C4H6 + C4H4"
  arrhenius_coeff = 1.4724e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1062 {
  chem_eq = "C6H5 + DIMEPTD --> 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7 +" &
        " 0.5*C4H6 + C6H6"
  arrhenius_coeff = 1.3320e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1063 {
  chem_eq = "C5H5 + DIMEPTD --> 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7 +" &
        " 0.5*C4H6 + C5H6"
  arrhenius_coeff = 1.2240e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1064 {
  chem_eq = "C3H3 + DIMEPTD --> 0.5*C3H5_A + 1.5*C3H4_A + 0.5*IC4H7 +" &
        " 0.5*C4H6"
  arrhenius_coeff = 9.7200e+04 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1065 {
  chem_eq = "LC5H7 + DIMEPTD --> 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7 +" &
        " 0.5*C4H6 + LC5H8"
  arrhenius_coeff = 9.7200e+04 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1066 {
  chem_eq = "DIMEPTD + C7H7 --> 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7 +" &
        " 0.5*C4H6 + C7H8"
  arrhenius_coeff = 4.8600e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1067 {
  chem_eq = "DIMEPTD + CH3C6H4 --> 0.5*C3H5_A + 0.5*C3H4_A + 0.5*IC4H7 +" &
        " 0.5*C4H6 + C7H8"
  arrhenius_coeff = 1.3320e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1068 {
  chem_eq = "C2H2 + C5H6 --> C7H8"
  arrhenius_coeff = 9.0000e+11 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1069 {
  chem_eq = "CH3 + C5H5CH3 --> H2 + H + C7H8"
  arrhenius_coeff = 4.0000e+11 0.000 7600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1070 {
  chem_eq = "C4H71_3 + C4H6 --> H2 + CH3 + C7H8"
  arrhenius_coeff = 6.0000e+11 0.000 13000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1071 {
  chem_eq = "2.0*IC4H7 --> H2 + H + CH3 + C7H8"
  arrhenius_coeff = 1.0000e+10 0.000 6000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1072 {
  chem_eq = "2.0*C4H71_3 --> H2 + H + CH3 + C7H8"
  arrhenius_coeff = 7.5000e+11 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1073 {
  chem_eq = "C3H3 + C4H6 --> H + C7H8"
  arrhenius_coeff = 2.0000e+11 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1073_rev {
  chem_eq = " H + C7H8 --> C3H3 + C4H6 "
  arrhenius_coeff = 2.0000e+11 0.000 5000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1074 {
  chem_eq = "C2H2 + C5H5 --> C7H7"
  arrhenius_coeff = 2.6100e+34 -13.610 37194.67
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 2.610000e+34" &
        " -1.361000e+01 3.719467e+04, 1.920000e+11 -6.850000e+00" &
        " 2.413656e+04, 1.050000e-02 -2.870000e+00 1.869178e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1074_rev {
  chem_eq = " C7H7 --> C2H2 + C5H5 "
  arrhenius_coeff = 2.6100e+34 -13.610 37194.67
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 arrhenius_press: 2.610000e+34" &
        " -1.361000e+01 3.719467e+04, 1.920000e+11 -6.850000e+00" &
        " 2.413656e+04, 1.050000e-02 -2.870000e+00 1.869178e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1075 {
  chem_eq = "C7H8 --> H + C7H7"
  arrhenius_coeff = 5.6000e+15 0.170 91168.00
  press_rxn_param = Troe_falloff "arrhenius_press: 1.00e+98" &
        " -22.855 99882.0 press_coeff: 0.06547 15.11 1.000e+10" &
        " 7.596e+07"
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1075_rev {
  chem_eq = " H + C7H7 --> C7H8 "
  arrhenius_coeff = 5.6000e+15 0.170 91168.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 1.00e+98" &
        " -22.855 99882.0 press_coeff: 0.06547 15.11 1.000e+10" &
        " 7.596e+07"
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1076 {
  chem_eq = "H + CH3C6H4 --> C7H8"
  arrhenius_coeff = 4.9000e+13 0.150 0.00
  press_rxn_param = Troe_falloff "arrhenius_press: 3.30e+75" &
        " -16.300 7000.0 press_coeff: 1.000 0.1000 584.9 6113."
  third_body_param = M_all "H2:2.00 CH4:2.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1076_rev {
  chem_eq = " C7H8 --> H + CH3C6H4 "
  arrhenius_coeff = 4.9000e+13 0.150 0.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 3.30e+75" &
        " -16.300 7000.0 press_coeff: 1.000 0.1000 584.9 6113."
  third_body_param = M_all "H2:2.00 CH4:2.00 "
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1077 {
  chem_eq = "C7H8 --> CH3 + C6H5"
  arrhenius_coeff = 1.9500e+27 -3.160 107447.00
  press_rxn_param = Troe_falloff "arrhenius_press: 1.00e+98" &
        " -22.966 122080.0 press_coeff: 0.7055 1.000e+10 459.9" &
        " 8.214e+09"
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1077_rev {
  chem_eq = " CH3 + C6H5 --> C7H8 "
  arrhenius_coeff = 1.9500e+27 -3.160 107447.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = Troe_falloff "arrhenius_press: 1.00e+98" &
        " -22.966 122080.0 press_coeff: 0.7055 1.000e+10 459.9" &
        " 8.214e+09"
  third_body_param = M_all
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1078 {
  chem_eq = "H + C7H7 --> CH3 + C6H5"
  arrhenius_coeff = 4.5000e+58 -11.900 51860.00
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.320000e-01 1.000000e+00 1.000000e+01 arrhenius_press:" &
        " 4.500000e+58 -1.190000e+01 5.186000e+04, 2.030000e+64" &
        " -1.337000e+01 5.952000e+04, 5.830000e+67 -1.415000e+01" &
        " 6.833000e+04, 8.850000e+68 -1.423000e+01 7.841000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1078_rev {
  chem_eq = " CH3 + C6H5 --> H + C7H7 "
  arrhenius_coeff = 4.5000e+58 -11.900 51860.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 3.950000e-02" &
        " 1.320000e-01 1.000000e+00 1.000000e+01 arrhenius_press:" &
        " 4.500000e+58 -1.190000e+01 5.186000e+04, 2.030000e+64" &
        " -1.337000e+01 5.952000e+04, 5.830000e+67 -1.415000e+01" &
        " 6.833000e+04, 8.850000e+68 -1.423000e+01 7.841000e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1079 {
  chem_eq = "CH3C6H4 --> C7H7"
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 5.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.050000e+65 -1.573000e+01 6.543690e+04," &
        " 5.580000e+45 -9.832000e+00 5.761785e+04, 7.600000e+24" &
        " -3.522000e+00 4.912007e+04, 2.790000e+14 -4.220000e-01" &
        " 4.368101e+04, 3.340000e+13 0.000000e+00 4.592354e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1079_rev {
  chem_eq = " C7H7 --> CH3C6H4 "
  arrhenius_coeff = 1.0000e+00 1.000 1.00
  reverse_calc = fromForwardRateConstant
  press_rxn_param = PLOG "press_coeff: 1.000000e-01" &
        " 1.000000e+00 1.000000e+01 5.000000e+01 1.000000e+02" &
        " arrhenius_press: 1.050000e+65 -1.573000e+01 6.543690e+04," &
        " 5.580000e+45 -9.832000e+00 5.761785e+04, 7.600000e+24" &
        " -3.522000e+00 4.912007e+04, 2.790000e+14 -4.220000e-01" &
        " 4.368101e+04, 3.340000e+13 0.000000e+00 4.592354e+04"
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1080 {
  chem_eq = "H + C7H8 --> CH3 + C6H6"
  arrhenius_coeff = 5.2800e+07 1.710 4387.13
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1080_rev {
  chem_eq = " CH3 + C6H6 --> H + C7H8 "
  arrhenius_coeff = 5.2800e+07 1.710 4387.13
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1081 {
  chem_eq = "H + C7H8 --> H2 + CH3C6H4"
  arrhenius_coeff = 4.1500e+07 2.120 14821.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1081_rev {
  chem_eq = " H2 + CH3C6H4 --> H + C7H8 "
  arrhenius_coeff = 4.1500e+07 2.120 14821.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1082 {
  chem_eq = "H + C7H8 --> H2 + C7H7"
  arrhenius_coeff = 3.1000e+04 2.760 3992.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1082_rev {
  chem_eq = " H2 + C7H7 --> H + C7H8 "
  arrhenius_coeff = 3.1000e+04 2.760 3992.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1083 {
  chem_eq = "C3H5_A + C4H2 --> C7H7"
  arrhenius_coeff = 2.0000e+11 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1084 {
  chem_eq = "C2H + C5H6 --> C7H7"
  arrhenius_coeff = 3.0000e+11 0.000 8000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1085 {
  chem_eq = "C7H7 --> C3H3 + C4H4"
  arrhenius_coeff = 1.0000e+14 0.000 83600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1086 {
  chem_eq = "C3H3 + C4H4 --> C7H7"
  arrhenius_coeff = 2.0000e+13 0.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1087 {
  chem_eq = "C3H4_P + C4H3 --> CH3C6H4"
  arrhenius_coeff = 2.4000e+10 0.000 0.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1088 {
  chem_eq = "CH3C6H4 --> C3H4_P + C4H3"
  arrhenius_coeff = 4.5000e+13 0.000 72500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1089 {
  chem_eq = "C2H4 + C4H6 --> H2 + CYC6H8"
  arrhenius_coeff = 2.0000e+12 0.000 36000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1090 {
  chem_eq = "C2H2 + C4H6 --> CYC6H8"
  arrhenius_coeff = 1.5000e+11 0.000 22500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1091 {
  chem_eq = "C3H6 + C4H6 --> H + CH3 + CYC6H8"
  arrhenius_coeff = 2.5000e+10 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1092 {
  chem_eq = "C3H5_A + IC4H7 --> H + CH3 + CYC6H8"
  arrhenius_coeff = 4.0000e+10 0.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1093 {
  chem_eq = "C3H5_S + C4H6 --> CH3 + CYC6H8"
  arrhenius_coeff = 2.0000e+11 0.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1094 {
  chem_eq = "C2H4 + C4H5 --> H + CYC6H8"
  arrhenius_coeff = 3.0000e+11 0.000 7000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1095 {
  chem_eq = "C4H6 + C4H5 --> C2H3 + CYC6H8"
  arrhenius_coeff = 5.0000e+11 0.000 7000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1096 {
  chem_eq = "H + CYC6H8 --> H2 + H + C6H6"
  arrhenius_coeff = 2.7000e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1097 {
  chem_eq = "CH3 + CYC6H8 --> H + CH4 + C6H6"
  arrhenius_coeff = 2.7000e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1098 {
  chem_eq = "C2H3 + CYC6H8 --> H + C2H4 + C6H6"
  arrhenius_coeff = 6.1020e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1099 {
  chem_eq = "C2H5 + CYC6H8 --> H + C2H6 + C6H6"
  arrhenius_coeff = 1.8000e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1100 {
  chem_eq = "NC3H7 + CYC6H8 --> H + C3H8 + C6H6"
  arrhenius_coeff = 1.2150e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1101 {
  chem_eq = "IC3H7 + CYC6H8 --> H + C3H8 + C6H6"
  arrhenius_coeff = 1.2150e+05 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1102 {
  chem_eq = "C3H5_S + CYC6H8 --> H + C3H6 + C6H6"
  arrhenius_coeff = 2.4300e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1103 {
  chem_eq = "C3H5_T + CYC6H8 --> H + C3H6 + C6H6"
  arrhenius_coeff = 2.4300e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1104 {
  chem_eq = "C3H5_A + CYC6H8 --> H + C3H6 + C6H6"
  arrhenius_coeff = 4.8510e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1105 {
  chem_eq = "C4H5 + CYC6H8 --> H + C4H6 + C6H6"
  arrhenius_coeff = 2.4300e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1106 {
  chem_eq = "C4H71_4 + CYC6H8 --> H + C4H8_1 + C6H6"
  arrhenius_coeff = 1.2150e+05 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1107 {
  chem_eq = "C4H71_3 + CYC6H8 --> H + C4H8_2 + C6H6"
  arrhenius_coeff = 3.0600e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1108 {
  chem_eq = "IC4H7 + CYC6H8 --> H + IC4H8 + C6H6"
  arrhenius_coeff = 2.4300e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1109 {
  chem_eq = "C4H3 + CYC6H8 --> H + C4H4 + C6H6"
  arrhenius_coeff = 3.6810e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1110 {
  chem_eq = "C6H5 + CYC6H8 --> H + 2.0*C6H6"
  arrhenius_coeff = 3.3300e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1111 {
  chem_eq = "C5H5 + CYC6H8 --> H + C6H6 + C5H6"
  arrhenius_coeff = 3.0600e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1112 {
  chem_eq = "C3H3 + CYC6H8 --> H + C3H4_A + C6H6"
  arrhenius_coeff = 2.4300e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1113 {
  chem_eq = "LC5H7 + CYC6H8 --> H + C6H6 + LC5H8"
  arrhenius_coeff = 2.4300e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1114 {
  chem_eq = "CYC6H8 + C7H7 --> H + C6H6 + C7H8"
  arrhenius_coeff = 1.2150e+05 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1115 {
  chem_eq = "CYC6H8 + CH3C6H4 --> H + C6H6 + C7H8"
  arrhenius_coeff = 3.3300e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1116 {
  chem_eq = "H + CYC6H10 --> H2 + 0.2*CH3 + 0.2*C5H6 + 0.8*RCYC6H9"
  arrhenius_coeff = 3.6000e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1117 {
  chem_eq = "CH3 + CYC6H10 --> CH4 + 0.2*CH3 + 0.2*C5H6 + 0.8*RCYC6H9"
  arrhenius_coeff = 3.6000e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1118 {
  chem_eq = "C2H3 + CYC6H10 --> 0.2*CH3 + C2H4 + 0.2*C5H6 + 0.8*RCYC6H9"
  arrhenius_coeff = 8.1360e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1119 {
  chem_eq = "C2H5 + CYC6H10 --> 0.2*CH3 + C2H6 + 0.2*C5H6 + 0.8*RCYC6H9"
  arrhenius_coeff = 2.4000e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1120 {
  chem_eq = "NC3H7 + CYC6H10 --> 0.2*CH3 + C3H8 + 0.2*C5H6 +" &
        " 0.8*RCYC6H9"
  arrhenius_coeff = 1.6200e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1121 {
  chem_eq = "IC3H7 + CYC6H10 --> 0.2*CH3 + C3H8 + 0.2*C5H6 +" &
        " 0.8*RCYC6H9"
  arrhenius_coeff = 1.6200e+05 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1122 {
  chem_eq = "C3H5_S + CYC6H10 --> 0.2*CH3 + C3H6 + 0.2*C5H6 +" &
        " 0.8*RCYC6H9"
  arrhenius_coeff = 3.2400e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1123 {
  chem_eq = "C3H5_T + CYC6H10 --> 0.2*CH3 + C3H6 + 0.2*C5H6 +" &
        " 0.8*RCYC6H9"
  arrhenius_coeff = 3.2400e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1124 {
  chem_eq = "C3H5_A + CYC6H10 --> 0.2*CH3 + C3H6 + 0.2*C5H6 +" &
        " 0.8*RCYC6H9"
  arrhenius_coeff = 6.4680e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1125 {
  chem_eq = "C4H5 + CYC6H10 --> 0.2*CH3 + C4H6 + 0.2*C5H6 + 0.8*RCYC6H9"
  arrhenius_coeff = 3.2400e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1126 {
  chem_eq = "C4H71_4 + CYC6H10 --> 0.2*CH3 + C4H8_1 + 0.2*C5H6 +" &
        " 0.8*RCYC6H9"
  arrhenius_coeff = 1.6200e+05 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1127 {
  chem_eq = "C4H71_3 + CYC6H10 --> 0.2*CH3 + C4H8_2 + 0.2*C5H6 +" &
        " 0.8*RCYC6H9"
  arrhenius_coeff = 4.0800e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1128 {
  chem_eq = "IC4H7 + CYC6H10 --> 0.2*CH3 + IC4H8 + 0.2*C5H6 +" &
        " 0.8*RCYC6H9"
  arrhenius_coeff = 3.2400e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1129 {
  chem_eq = "C4H3 + CYC6H10 --> 0.2*CH3 + C4H4 + 0.2*C5H6 + 0.8*RCYC6H9"
  arrhenius_coeff = 4.9080e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1130 {
  chem_eq = "C6H5 + CYC6H10 --> 0.2*CH3 + C6H6 + 0.2*C5H6 + 0.8*RCYC6H9"
  arrhenius_coeff = 4.4400e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1131 {
  chem_eq = "C5H5 + CYC6H10 --> 0.2*CH3 + 1.2*C5H6 + 0.8*RCYC6H9"
  arrhenius_coeff = 4.0800e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1132 {
  chem_eq = "C3H3 + CYC6H10 --> 0.2*CH3 + C3H4_A + 0.2*C5H6 +" &
        " 0.8*RCYC6H9"
  arrhenius_coeff = 3.2400e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1133 {
  chem_eq = "LC5H7 + CYC6H10 --> 0.2*CH3 + 0.2*C5H6 + LC5H8 +" &
        " 0.8*RCYC6H9"
  arrhenius_coeff = 3.2400e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1134 {
  chem_eq = "CYC6H10 + C7H7 --> 0.2*CH3 + 0.2*C5H6 + 0.8*RCYC6H9 + C7H8"
  arrhenius_coeff = 1.6200e+05 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1135 {
  chem_eq = "CYC6H10 + CH3C6H4 --> 0.2*CH3 + 0.2*C5H6 + 0.8*RCYC6H9 + C7H8"
  arrhenius_coeff = 4.4400e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1136 {
  chem_eq = "CH3 + C7H8 --> CH4 + C7H7"
  arrhenius_coeff = 9.0000e+04 2.000 5948.43
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1137 {
  chem_eq = "C2H3 + C7H8 --> C2H4 + C7H7"
  arrhenius_coeff = 2.0340e+05 2.000 4890.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1138 {
  chem_eq = "C2H5 + C7H8 --> C2H6 + C7H7"
  arrhenius_coeff = 6.0000e+04 2.000 7090.53
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1139 {
  chem_eq = "NC3H7 + C7H8 --> C3H8 + C7H7"
  arrhenius_coeff = 4.0500e+04 2.000 7090.53
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1140 {
  chem_eq = "IC3H7 + C7H8 --> C3H8 + C7H7"
  arrhenius_coeff = 4.0500e+04 2.000 10290.53
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1141 {
  chem_eq = "C3H5_S + C7H8 --> C3H6 + C7H7"
  arrhenius_coeff = 8.1000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1142 {
  chem_eq = "C3H5_T + C7H8 --> C3H6 + C7H7"
  arrhenius_coeff = 8.1000e+04 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1143 {
  chem_eq = "C3H5_A + C7H8 --> C3H6 + C7H7"
  arrhenius_coeff = 1.6170e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1144 {
  chem_eq = "C4H5 + C7H8 --> C4H6 + C7H7"
  arrhenius_coeff = 8.1000e+04 2.000 5490.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1145 {
  chem_eq = "C4H71_4 + C7H8 --> C4H8_1 + C7H7"
  arrhenius_coeff = 4.0500e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1146 {
  chem_eq = "C4H71_3 + C7H8 --> C4H8_2 + C7H7"
  arrhenius_coeff = 1.0200e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1147 {
  chem_eq = "IC4H7 + C7H8 --> IC4H8 + C7H7"
  arrhenius_coeff = 8.1000e+04 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1148 {
  chem_eq = "C4H3 + C7H8 --> C4H4 + C7H7"
  arrhenius_coeff = 1.2270e+05 2.000 6148.43
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1149 {
  chem_eq = "C6H5 + C7H8 --> C6H6 + C7H7"
  arrhenius_coeff = 1.1100e+08 1.000 2990.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1150 {
  chem_eq = "C5H5 + C7H8 --> C5H6 + C7H7"
  arrhenius_coeff = 1.0200e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1151 {
  chem_eq = "C3H3 + C7H8 --> C3H4_A + C7H7"
  arrhenius_coeff = 8.1000e+04 2.000 14990.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1152 {
  chem_eq = "LC5H7 + C7H8 --> LC5H8 + C7H7"
  arrhenius_coeff = 8.1000e+04 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1153 {
  chem_eq = "C7H8 + CH3C6H4 --> C7H8 + C7H7"
  arrhenius_coeff = 1.1100e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1154 {
  chem_eq = "CH3 + C7H8 --> CH4 + CH3C6H4"
  arrhenius_coeff = 3.0000e+05 2.000 8925.13
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1155 {
  chem_eq = "C2H3 + C7H8 --> C2H4 + CH3C6H4"
  arrhenius_coeff = 6.7800e+05 2.000 7460.65
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1156 {
  chem_eq = "C2H5 + C7H8 --> C2H6 + CH3C6H4"
  arrhenius_coeff = 2.0000e+05 2.000 12330.09
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1157 {
  chem_eq = "NC3H7 + C7H8 --> C3H8 + CH3C6H4"
  arrhenius_coeff = 1.3500e+05 2.000 12330.09
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1158 {
  chem_eq = "IC3H7 + C7H8 --> C3H8 + CH3C6H4"
  arrhenius_coeff = 1.3500e+05 2.000 16129.68
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1159 {
  chem_eq = "C3H5_S + C7H8 --> C3H6 + CH3C6H4"
  arrhenius_coeff = 2.7000e+05 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1160 {
  chem_eq = "C3H5_T + C7H8 --> C3H6 + CH3C6H4"
  arrhenius_coeff = 2.7000e+05 2.000 10579.25
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1161 {
  chem_eq = "C3H5_A + C7H8 --> C3H6 + CH3C6H4"
  arrhenius_coeff = 5.3900e+05 2.000 20432.90
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1162 {
  chem_eq = "C4H5 + C7H8 --> C4H6 + CH3C6H4"
  arrhenius_coeff = 2.7000e+05 2.000 8207.64
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1163 {
  chem_eq = "C4H71_4 + C7H8 --> C4H8_1 + CH3C6H4"
  arrhenius_coeff = 1.3500e+05 2.000 13144.81
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1164 {
  chem_eq = "C4H71_3 + C7H8 --> C4H8_2 + CH3C6H4"
  arrhenius_coeff = 3.4000e+05 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1165 {
  chem_eq = "IC4H7 + C7H8 --> IC4H8 + CH3C6H4"
  arrhenius_coeff = 2.7000e+05 2.000 21467.03
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1166 {
  chem_eq = "C4H3 + C7H8 --> C4H4 + CH3C6H4"
  arrhenius_coeff = 4.0900e+05 2.000 9170.82
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1167 {
  chem_eq = "C6H5 + C7H8 --> C6H6 + CH3C6H4"
  arrhenius_coeff = 3.7000e+08 1.000 4356.78
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1168 {
  chem_eq = "C5H5 + C7H8 --> C5H6 + CH3C6H4"
  arrhenius_coeff = 3.4000e+05 2.000 21810.72
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1169 {
  chem_eq = "C3H3 + C7H8 --> C3H4_A + CH3C6H4"
  arrhenius_coeff = 2.7000e+05 2.000 19439.11
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1170 {
  chem_eq = "LC5H7 + C7H8 --> LC5H8 + CH3C6H4"
  arrhenius_coeff = 2.7000e+05 2.000 22952.95
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1171 {
  chem_eq = "C7H8 + C7H7 --> C7H8 + CH3C6H4"
  arrhenius_coeff = 1.3500e+05 2.000 24090.27
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1172 {
  chem_eq = "C16H32 --> CH3 + C2H4 + C3H6 + C3H5_A + NC7H14"
  arrhenius_coeff = 1.5000e+16 0.000 71000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1173 {
  chem_eq = "C16H32 --> CH3 + C2H4 + C3H6 + C3H5_A + NC7H14"
  arrhenius_coeff = 3.0000e+18 0.000 86000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1174 {
  chem_eq = "C16H32 --> H + C2H5 + C3H6 + C4H6 + NC7H14"
  arrhenius_coeff = 3.0000e+18 0.000 86000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1175 {
  chem_eq = "H + C16H32 --> H2 + C2H5 + C3H6 + C4H6 + NC7H14"
  arrhenius_coeff = 1.9800e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1176 {
  chem_eq = "CH3 + C16H32 --> CH4 + C2H5 + C3H6 + C4H6 + NC7H14"
  arrhenius_coeff = 1.9800e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1177 {
  chem_eq = "C2H3 + C16H32 --> C2H5 + C2H4 + C3H6 + C4H6 + NC7H14"
  arrhenius_coeff = 4.4748e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1178 {
  chem_eq = "C2H5 + C16H32 --> C2H6 + C2H5 + C3H6 + C4H6 + NC7H14"
  arrhenius_coeff = 1.3200e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1179 {
  chem_eq = "NC3H7 + C16H32 --> C2H5 + C3H8 + C3H6 + C4H6 + NC7H14"
  arrhenius_coeff = 8.9100e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1180 {
  chem_eq = "IC3H7 + C16H32 --> C2H5 + C3H8 + C3H6 + C4H6 + NC7H14"
  arrhenius_coeff = 8.9100e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1181 {
  chem_eq = "C3H5_S + C16H32 --> C2H5 + 2.0*C3H6 + C4H6 + NC7H14"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1182 {
  chem_eq = "C3H5_T + C16H32 --> C2H5 + 2.0*C3H6 + C4H6 + NC7H14"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1183 {
  chem_eq = "C3H5_A + C16H32 --> C2H5 + 2.0*C3H6 + C4H6 + NC7H14"
  arrhenius_coeff = 3.5574e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1184 {
  chem_eq = "C4H5 + C16H32 --> C2H5 + C3H6 + 2.0*C4H6 + NC7H14"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1185 {
  chem_eq = "C4H71_4 + C16H32 --> C2H5 + C3H6 + C4H8_1 + C4H6 + NC7H14"
  arrhenius_coeff = 8.9100e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1186 {
  chem_eq = "C4H71_3 + C16H32 --> C2H5 + C3H6 + C4H8_2 + C4H6 + NC7H14"
  arrhenius_coeff = 2.2440e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1187 {
  chem_eq = "IC4H7 + C16H32 --> C2H5 + C3H6 + IC4H8 + C4H6 + NC7H14"
  arrhenius_coeff = 1.7820e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1188 {
  chem_eq = "C4H3 + C16H32 --> C2H5 + C3H6 + C4H6 + C4H4 + NC7H14"
  arrhenius_coeff = 2.6994e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1189 {
  chem_eq = "C6H5 + C16H32 --> C2H5 + C3H6 + C4H6 + C6H6 + NC7H14"
  arrhenius_coeff = 2.4420e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1190 {
  chem_eq = "C5H5 + C16H32 --> C2H5 + C3H6 + C4H6 + C5H6 + NC7H14"
  arrhenius_coeff = 2.2440e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1191 {
  chem_eq = "C3H3 + C16H32 --> C2H5 + C3H6 + C3H4_A + C4H6 + NC7H14"
  arrhenius_coeff = 1.7820e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1192 {
  chem_eq = "LC5H7 + C16H32 --> C2H5 + C3H6 + C4H6 + LC5H8 + NC7H14"
  arrhenius_coeff = 1.7820e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1193 {
  chem_eq = "C7H7 + C16H32 --> C2H5 + C3H6 + C4H6 + NC7H14 + C7H8"
  arrhenius_coeff = 8.9100e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1194 {
  chem_eq = "CH3C6H4 + C16H32 --> C2H5 + C3H6 + C4H6 + NC7H14 + C7H8"
  arrhenius_coeff = 2.4420e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1195 {
  chem_eq = "H + C16H32 --> H2 + 0.44*H + 0.08*CH3 + 1.52*C2H4 +" &
        " 0.8*C3H6 + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 7.5000e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1196 {
  chem_eq = "CH3 + C16H32 --> 0.44*H + CH4 + 0.08*CH3 + 1.52*C2H4 +" &
        " 0.8*C3H6 + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 7.5000e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1197 {
  chem_eq = "C2H3 + C16H32 --> 0.44*H + 0.08*CH3 + 2.52*C2H4 + 0.8*C3H6" &
        " + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 + 0.08*LC5H8 +" &
        " 0.8*NC7H14"
  arrhenius_coeff = 1.6950e+06 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1198 {
  chem_eq = "C2H5 + C16H32 --> 0.44*H + 0.08*CH3 + C2H6 + 1.52*C2H4 +" &
        " 0.8*C3H6 + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 5.0000e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1199 {
  chem_eq = "NC3H7 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 + C3H8 +" &
        " 0.8*C3H6 + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 3.3750e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1200 {
  chem_eq = "IC3H7 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 + C3H8 +" &
        " 0.8*C3H6 + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 3.3750e+05 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1201 {
  chem_eq = "C3H5_S + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 +" &
        " 1.8*C3H6 + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 6.7500e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1202 {
  chem_eq = "C3H5_T + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 +" &
        " 1.8*C3H6 + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 6.7500e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1203 {
  chem_eq = "C3H5_A + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 +" &
        " 1.8*C3H6 + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 1.3475e+06 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1204 {
  chem_eq = "C4H5 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 + 0.8*C3H6" &
        " + 0.48*C3H5_A + 0.32*C4H8_1 + 1.44*C4H6 + 0.08*LC5H8 +" &
        " 0.8*NC7H14"
  arrhenius_coeff = 6.7500e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1205 {
  chem_eq = "C4H71_4 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 +" &
        " 0.8*C3H6 + 0.48*C3H5_A + 1.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 3.3750e+05 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1206 {
  chem_eq = "C4H71_3 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 +" &
        " 0.8*C3H6 + 0.48*C3H5_A + 0.32*C4H8_1 + C4H8_2 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 8.5000e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1207 {
  chem_eq = "IC4H7 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 + 0.8*C3H6" &
        " + 0.48*C3H5_A + IC4H8 + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 6.7500e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1208 {
  chem_eq = "C4H3 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 + 0.8*C3H6" &
        " + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 + C4H4 + 0.08*LC5H8" &
        " + 0.8*NC7H14"
  arrhenius_coeff = 1.0225e+06 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1209 {
  chem_eq = "C6H5 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 + 0.8*C3H6" &
        " + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 + C6H6 + 0.08*LC5H8" &
        " + 0.8*NC7H14"
  arrhenius_coeff = 9.2500e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1210 {
  chem_eq = "C5H5 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 + 0.8*C3H6" &
        " + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 + C5H6 + 0.08*LC5H8" &
        " + 0.8*NC7H14"
  arrhenius_coeff = 8.5000e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1211 {
  chem_eq = "C3H3 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 + 0.8*C3H6" &
        " + 0.48*C3H5_A + C3H4_A + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14"
  arrhenius_coeff = 6.7500e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1212 {
  chem_eq = "LC5H7 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 + 0.8*C3H6" &
        " + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 + 1.08*LC5H8 +" &
        " 0.8*NC7H14"
  arrhenius_coeff = 6.7500e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1213 {
  chem_eq = "C7H7 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 + 0.8*C3H6" &
        " + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 + 0.08*LC5H8 +" &
        " 0.8*NC7H14 + C7H8"
  arrhenius_coeff = 3.3750e+05 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1214 {
  chem_eq = "CH3C6H4 + C16H32 --> 0.44*H + 0.08*CH3 + 1.52*C2H4 +" &
        " 0.8*C3H6 + 0.48*C3H5_A + 0.32*C4H8_1 + 0.44*C4H6 +" &
        " 0.08*LC5H8 + 0.8*NC7H14 + C7H8"
  arrhenius_coeff = 9.2500e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1215 {
  chem_eq = "C30H60 --> CH3 + 5.0*C2H4 + C3H5_A + C16H32"
  arrhenius_coeff = 1.5000e+16 0.000 71000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1216 {
  chem_eq = "C30H60 --> CH3 + 5.0*C2H4 + C3H5_A + C16H32"
  arrhenius_coeff = 6.0000e+18 0.000 86000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1217 {
  chem_eq = "C30H60 --> H + C2H5 + 4.0*C2H4 + C4H6 + C16H32"
  arrhenius_coeff = 6.0000e+18 0.000 86000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1218 {
  chem_eq = "H + C30H60 --> H2 + CH3 + C2H4 + C4H6 + NC7H14 + C16H32"
  arrhenius_coeff = 1.9800e+07 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1219 {
  chem_eq = "CH3 + C30H60 --> CH4 + CH3 + C2H4 + C4H6 + NC7H14 + C16H32"
  arrhenius_coeff = 1.9800e+05 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1220 {
  chem_eq = "C2H3 + C30H60 --> CH3 + 2.0*C2H4 + C4H6 + NC7H14 + C16H32"
  arrhenius_coeff = 4.4748e+05 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1221 {
  chem_eq = "C2H5 + C30H60 --> CH3 + C2H6 + C2H4 + C4H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 1.3200e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1222 {
  chem_eq = "NC3H7 + C30H60 --> CH3 + C2H4 + C3H8 + C4H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 8.9100e+04 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1223 {
  chem_eq = "IC3H7 + C30H60 --> CH3 + C2H4 + C3H8 + C4H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 8.9100e+04 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1224 {
  chem_eq = "C3H5_S + C30H60 --> CH3 + C2H4 + C3H6 + C4H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1225 {
  chem_eq = "C3H5_T + C30H60 --> CH3 + C2H4 + C3H6 + C4H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1226 {
  chem_eq = "C3H5_A + C30H60 --> CH3 + C2H4 + C3H6 + C4H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 3.5574e+05 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1227 {
  chem_eq = "C4H5 + C30H60 --> CH3 + C2H4 + 2.0*C4H6 + NC7H14 + C16H32"
  arrhenius_coeff = 1.7820e+05 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1228 {
  chem_eq = "C4H71_4 + C30H60 --> CH3 + C2H4 + C4H8_1 + C4H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 8.9100e+04 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1229 {
  chem_eq = "C4H71_3 + C30H60 --> CH3 + C2H4 + C4H8_2 + C4H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 2.2440e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1230 {
  chem_eq = "IC4H7 + C30H60 --> CH3 + C2H4 + IC4H8 + C4H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 1.7820e+05 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1231 {
  chem_eq = "C4H3 + C30H60 --> CH3 + C2H4 + C4H6 + C4H4 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 2.6994e+05 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1232 {
  chem_eq = "C6H5 + C30H60 --> CH3 + C2H4 + C4H6 + C6H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 2.4420e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1233 {
  chem_eq = "C5H5 + C30H60 --> CH3 + C2H4 + C4H6 + C5H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 2.2440e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1234 {
  chem_eq = "C3H3 + C30H60 --> CH3 + C2H4 + C3H4_A + C4H6 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 1.7820e+05 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1235 {
  chem_eq = "LC5H7 + C30H60 --> CH3 + C2H4 + C4H6 + LC5H8 + NC7H14 +" &
        " C16H32"
  arrhenius_coeff = 1.7820e+05 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1236 {
  chem_eq = "C7H7 + C30H60 --> CH3 + C2H4 + C4H6 + NC7H14 + C7H8 +" &
        " C16H32"
  arrhenius_coeff = 8.9100e+04 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1237 {
  chem_eq = "CH3C6H4 + C30H60 --> CH3 + C2H4 + C4H6 + NC7H14 + C7H8 +" &
        " C16H32"
  arrhenius_coeff = 2.4420e+08 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1238 {
  chem_eq = "H + C30H60 --> 0.999955348837208*H2 + 0.469977674418604*H +" &
        " 0.039960930232558*CH3 + 1.641966511627907*C2H4 +" &
        " 0.48994976744186*C3H6 + 0.489972093023256*C3H5_A +" &
        " 0.259933023255814*C4H8_1 + 0.469977674418606*C4H6 +" &
        " 0.039960930232559*LC5H8 + 0.529882790697674*NC7H14 +" &
        " 1.056732093023256*C16H32"
  arrhenius_coeff = 1.5900e+08 2.000 4000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1239 {
  chem_eq = "CH3 + C30H60 --> 0.469980241492864*H +" &
        " 0.999944017563116*CH4 + 0.039963776070252*CH3 +" &
        " 1.641967069154775*C2H4 + 0.489950603732162*C3H6 +" &
        " 0.489970362239298*C3H5_A + 0.25993413830955*C4H8_1 +" &
        " 0.469973655323821*C4H6 + 0.039957189901208*LC5H8 +" &
        " 0.529884742041712*NC7H14 + 1.0567365532382*C16H32"
  arrhenius_coeff = 1.5900e+06 2.000 5000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1240 {
  chem_eq = "C2H3 + C30H60 --> 0.469968709256845*H +" &
        " 0.03995149934811*CH3 + 2.641965580182529*C2H4 +" &
        " 0.489948370273794*C3H6 + 0.489979661016949*C3H5_A +" &
        " 0.259931160365059*C4H8_1 + 0.469993741851369*C4H6 +" &
        " 0.039976531942634*LC5H8 + 0.529879530638853*NC7H14 +" &
        " 1.056724641460235*C16H32"
  arrhenius_coeff = 3.5934e+06 2.000 3900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1241 {
  chem_eq = "C2H5 + C30H60 --> 0.469982792615164*H +" &
        " 0.03996666069188*CH3 + 0.999933321383761*C2H6 +" &
        " 1.641967736153432*C2H4 + 0.489951604230149*C3H6 +" &
        " 0.489968811614985*C3H5_A + 0.259935472306865*C4H8_1 +" &
        " 0.469969887076537*C4H6 + 0.039953755153253*LC5H8 +" &
        " 0.529887076537014*NC7H14 + 1.05674188922746*C16H32"
  arrhenius_coeff = 1.0600e+06 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1242 {
  chem_eq = "NC3H7 + C30H60 --> 0.469985299264963*H +" &
        " 0.039969548477424*CH3 + 1.641968498424921*C2H4 +" &
        " 0.999923346167308*C3H8 + 0.489952747637382*C3H6 +" &
        " 0.489967448372419*C3H5_A + 0.259936996849843*C4H8_1 +" &
        " 0.469966398319916*C4H6 + 0.039950647532377*LC5H8 +" &
        " 0.529889744487224*NC7H14 + 1.05674798739937*C16H32"
  arrhenius_coeff = 7.1550e+05 2.000 6700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1243 {
  chem_eq = "IC3H7 + C30H60 --> 0.469985299264963*H +" &
        " 0.039969548477424*CH3 + 1.641968498424921*C2H4 +" &
        " 0.999923346167308*C3H8 + 0.489952747637382*C3H6 +" &
        " 0.489967448372419*C3H5_A + 0.259936996849843*C4H8_1 +" &
        " 0.469966398319916*C4H6 + 0.039950647532377*LC5H8 +" &
        " 0.529889744487224*NC7H14 + 1.05674798739937*C16H32"
  arrhenius_coeff = 7.1550e+05 2.000 9900.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1244 {
  chem_eq = "C3H5_S + C30H60 --> 0.469968709256845*H +" &
        " 0.03995149934811*CH3 + 1.641965580182529*C2H4 +" &
        " 1.489948370273794*C3H6 + 0.489979661016949*C3H5_A +" &
        " 0.259931160365059*C4H8_1 + 0.469993741851369*C4H6 +" &
        " 0.039976531942634*LC5H8 + 0.529879530638853*NC7H14 +" &
        " 1.056724641460235*C16H32"
  arrhenius_coeff = 1.4310e+06 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1245 {
  chem_eq = "C3H5_T + C30H60 --> 0.469968709256845*H +" &
        " 0.03995149934811*CH3 + 1.641965580182529*C2H4 +" &
        " 1.489948370273794*C3H6 + 0.489979661016949*C3H5_A +" &
        " 0.259931160365059*C4H8_1 + 0.469993741851369*C4H6 +" &
        " 0.039976531942634*LC5H8 + 0.529879530638853*NC7H14 +" &
        " 1.056724641460235*C16H32"
  arrhenius_coeff = 1.4310e+06 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1246 {
  chem_eq = "C3H5_A + C30H60 --> 0.469968709256845*H +" &
        " 0.03995149934811*CH3 + 1.641965580182529*C2H4 +" &
        " 1.489948370273794*C3H6 + 0.489979661016949*C3H5_A +" &
        " 0.259931160365059*C4H8_1 + 0.469993741851369*C4H6 +" &
        " 0.039976531942634*LC5H8 + 0.529879530638853*NC7H14 +" &
        " 1.056724641460235*C16H32"
  arrhenius_coeff = 2.8567e+06 2.000 12800.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1247 {
  chem_eq = "C4H5 + C30H60 --> 0.469968709256845*H +" &
        " 0.03995149934811*CH3 + 1.641965580182529*C2H4 +" &
        " 0.489948370273794*C3H6 + 0.489979661016949*C3H5_A +" &
        " 0.259931160365059*C4H8_1 + 1.469993741851369*C4H6 +" &
        " 0.039976531942634*LC5H8 + 0.529879530638853*NC7H14 +" &
        " 1.056724641460235*C16H32"
  arrhenius_coeff = 1.4310e+06 2.000 4500.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1248 {
  chem_eq = "C4H71_4 + C30H60 --> 0.469968709256845*H +" &
        " 0.03995149934811*CH3 + 1.641965580182529*C2H4 +" &
        " 0.489948370273794*C3H6 + 0.489979661016949*C3H5_A +" &
        " 1.259931160365059*C4H8_1 + 0.469993741851369*C4H6 +" &
        " 0.039976531942634*LC5H8 + 0.529879530638853*NC7H14 +" &
        " 1.056724641460235*C16H32"
  arrhenius_coeff = 7.1550e+05 2.000 6600.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1249 {
  chem_eq = "C4H71_3 + C30H60 --> 0.469970082273747*H +" &
        " 0.039953627524308*CH3 + 1.641967090501122*C2H4 +" &
        " 0.489950635751683*C3H6 + 0.489980553477936*C3H5_A +" &
        " 0.259934181002244*C4H8_1 + 0.999934181002244*C4H8_2 +" &
        " 0.469994016454749*C4H6 + 0.03997756170531*LC5H8 +" &
        " 0.529884816753927*NC7H14 + 1.056736724008975*C16H32"
  arrhenius_coeff = 1.8020e+06 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1250 {
  chem_eq = "IC4H7 + C30H60 --> 0.469970082273747*H +" &
        " 0.039953627524308*CH3 + 1.641967090501122*C2H4 +" &
        " 0.489950635751683*C3H6 + 0.489980553477936*C3H5_A +" &
        " 0.999934181002244*IC4H8 + 0.259934181002244*C4H8_1 +" &
        " 0.469994016454749*C4H6 + 0.03997756170531*LC5H8 +" &
        " 0.529884816753927*NC7H14 + 1.056736724008975*C16H32"
  arrhenius_coeff = 1.4310e+06 2.000 13700.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1251 {
  chem_eq = "C4H3 + C30H60 --> 0.469977337110482*H +" &
        " 0.039960339943343*CH3 + 1.641966005665722*C2H4 +" &
        " 0.489949008498584*C3H6 + 0.489971671388102*C3H5_A +" &
        " 0.259932011331445*C4H8_1 + 0.469977337110481*C4H6 +" &
        " 1.000022662889517*C4H4 + 0.039960339943342*LC5H8 +" &
        " 0.529881019830029*NC7H14 + 1.056728045325779*C16H32"
  arrhenius_coeff = 2.1677e+06 2.000 5200.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1252 {
  chem_eq = "C6H5 + C30H60 --> 0.46997982103177*H +" &
        " 0.039962885112005*CH3 + 1.641966128160471*C2H4 +" &
        " 0.489949192240706*C3H6 + 0.489969371208936*C3H5_A +" &
        " 0.259932256320942*C4H8_1 + 0.469972614257402*C4H6 +" &
        " 1.000019458290793*C6H6 + 0.039955678337637*LC5H8 +" &
        " 0.529881448561648*NC7H14 + 1.056729025283767*C16H32"
  arrhenius_coeff = 1.9610e+09 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1253 {
  chem_eq = "C5H5 + C30H60 --> 0.469974656810982*H +" &
        " 0.039957550158395*CH3 + 1.641965786694826*C2H4 +" &
        " 0.489948680042239*C3H6 + 0.489974023231257*C3H5_A +" &
        " 0.259931573389652*C4H8_1 + 0.469982259767687*C4H6 +" &
        " 1.000015839493136*C5H6 + 0.0399651531151*LC5H8 +" &
        " 0.52988025343189*NC7H14 + 1.056726293558606*C16H32"
  arrhenius_coeff = 1.8020e+06 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1254 {
  chem_eq = "C3H3 + C30H60 --> 0.469970194879633*H +" &
        " 0.039952999617883*CH3 + 1.6419656094765*C2H4 +" &
        " 0.48994841421475*C3H6 + 0.489978219335117*C3H5_A +" &
        " 1.000008024455483*C3H4_A + 0.259931218953*C4H8_1 +" &
        " 0.469990829193733*C4H6 + 0.039973633931983*LC5H8 +" &
        " 0.529879633167749*NC7H14 + 1.056724875811998*C16H32"
  arrhenius_coeff = 1.4310e+06 2.000 14000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1255 {
  chem_eq = "LC5H7 + C30H60 --> 0.469968709256845*H +" &
        " 0.03995149934811*CH3 + 1.641965580182529*C2H4 +" &
        " 0.489948370273794*C3H6 + 0.489979661016949*C3H5_A +" &
        " 0.259931160365059*C4H8_1 + 0.469993741851369*C4H6 +" &
        " 1.039976531942634*LC5H8 + 0.529879530638853*NC7H14 +" &
        " 1.056724641460235*C16H32"
  arrhenius_coeff = 1.4310e+06 2.000 15000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1256 {
  chem_eq = "C7H7 + C30H60 --> 0.469977528089888*H +" &
        " 0.039960492932222*CH3 + 1.641965929684668*C2H4 +" &
        " 0.489948894527003*C3H6 + 0.489971366437115*C3H5_A +" &
        " 0.259931859369337*C4H8_1 + 0.469976803189561*C4H6 +" &
        " 0.039959768031896*LC5H8 + 0.529880753896339*NC7H14 +" &
        " 1.000015585357013*C7H8 + 1.056727437477347*C16H32"
  arrhenius_coeff = 7.1550e+05 2.000 16000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1257 {
  chem_eq = "CH3C6H4 + C30H60 --> 0.469977528089888*H +" &
        " 0.039960492932222*CH3 + 1.641965929684668*C2H4 +" &
        " 0.489948894527003*C3H6 + 0.489971366437115*C3H5_A +" &
        " 0.259931859369337*C4H8_1 + 0.469976803189561*C4H6 +" &
        " 0.039959768031896*LC5H8 + 0.529880753896339*NC7H14 +" &
        " 1.000015585357013*C7H8 + 1.056727437477347*C16H32"
  arrhenius_coeff = 1.9610e+09 1.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1258 {
  chem_eq = "C3H4_A + C3H3 --> H + C6H6"
  arrhenius_coeff = 2.2000e+11 0.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1258_rev {
  chem_eq = " H + C6H6 --> C3H4_A + C3H3 "
  arrhenius_coeff = 2.2000e+11 0.000 2000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1259 {
  chem_eq = "C3H4_P + C3H3 --> H + C6H6"
  arrhenius_coeff = 2.2000e+11 0.000 2000.00
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1259_rev {
  chem_eq = " H + C6H6 --> C3H4_P + C3H3 "
  arrhenius_coeff = 2.2000e+11 0.000 2000.00
  reverse_calc = fromForwardRateConstant
  dh = 0
  fracdh(0) = 1.0
 }

Reaction_1260 {
  chem_eq = "C2H4 + C4H2 --> C6H6"
  arrhenius_coeff = 3.0000e+11 0.000 30000.00
  dh = 0
  fracdh(0) = 1.0
 }


"""
