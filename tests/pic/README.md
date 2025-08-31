# Test cases for MFIX PIC

## Test case overview

| Test  | Description                                        |
| ----  | -------------------------------------------------- |
| PIC01 | Terminal velocity                            	     |
| PIC02 | Advection - velocity interpolation                 |
| PIC03 | Advection - parcel volume deposition               |
| PIC04 | Particle-settling in a fluid                       |
| PIC05 | Evaporation                         		     |
| PIC06 | Rayleigh-Taylor instability                        |
| PIC07 | Minimum fluidization				     |



## Feature/model test overview
|              | PIC01 | PIC02 | PIC03 | PIC04 | PIC05 | PIC06 | PIC07 |
| ------------ | :---: | :---: | :---: | :---: | :---: | :---: | :---: | 
| Test Freq.   |   C   |   C   |   C   |   C   |   C   |   C   |   C   | 
| Ref data set |   A   |   A   |   A   |   A   |   P   |   P   |   P   | 
| Dimension    |  3D   |  3D   |  3D   |  3D   |  3D   |  3D   |  3D   |
| Coupled      |   X   |   X   |   X   |   X   |       |   X   |   X   | 
| Momentum     |       |   X   |   X   |   X   |       |   X   |   X   | 
| Them Energy  |       |       |       |       |   X   |       |       | 
| Species Mass |       |       |       |       |   X   |       |       | 
| Turbulence   |       |       |       |       |       |       |       | 
| SMP Parallel |       |       |       |       |       |       |       | 
| DMP Parallel |       |  X    |   X   |   X   |       |   X   |   X   | 



--------------------------------------------------------------------

## Notice
Neither the United States Government nor any agency thereof, nor any
of their employees, makes any warranty, expressed or implied, or
assumes any legal liability or responsibility for the accuracy,
completeness, or usefulness of any information, apparatus, product,
or process disclosed or represents that its use would not infringe
privately owned rights.

* MFIX is provided without any user support for applications in the
  user's immediate organization. It should not be redistributed in
  whole or in part.

* The use of MFIX is to be acknowledged in any published paper based
  on computations using this software by citing the MFIX theory
  manual. Some of the submodels are being developed by researchers
  outside of NETL. The use of such submodels is to be acknowledged
  by citing the appropriate papers of the developers of the submodels.

* The authors would appreciate receiving any reports of bugs or other
  difficulties with the software, enhancements to the software, and
  accounts of practical applications of this software.

# Disclaimer
This report was prepared as an account of work sponsored by an agency
of the United States Government. Neither the United States Government
nor any agency thereof, nor any of their employees, makes any
warranty, express or implied, or assumes any legal liability or
responsibility for the accuracy, completeness, or usefulness of any
information, apparatus, product, or process disclosed, or represents
that its use would not infringe privately owned rights. Reference
herein to any specific commercial product, process, or service by
trade name, trademark, manufacturer, or otherwise does not
necessarily constitute or imply its endorsement, recommendation, or
favoring by the United States Government or any agency thereof. The
views and opinions of authors expressed herein do not necessarily
state or reflect those of the United States Government or any
agency thereof.
