gfortran -c data.f90
gfortran -c functions.f90
gfortran -c overlap.f90

gfortran -c onsager_B2.f90
gfortran -c onsager_B3.f90
gfortran -c onsager_B4.f90
gfortran -c onsager_B5.f90

python3.5 -m numpy.f2py -c -m virial virial.f90 data.f90 functions.f90 overlap.f90 onsager_B2.f90 onsager_B3.f90 onsager_B4.f90 onsager_B5.f90
