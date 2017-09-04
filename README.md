# Spherocylinders
Virial coefficients for hard spherocylinders
 
Files:
- Bn_L5d0_D1d0-norm.dat     -> nth virial data in tab deliminated format
- Bn_L5d0_D1d0-norm-csv.dat -> nth virial data in comma separated format
- process_virials.ipynb     -> Notebook to analyse the data
- data_B5_Cc.dat            -> Contains the EOS data for C=c, for both the isotropic and nematic phases

The virial files contain:
1) index
2) S_nem: nematic order parameter
3) alpha (for S_nem = 1, alpha is set to be some very high value, since it is infinity for perfect alignment)
4) Bn: the value of the virial
5) error: the statistical error in the MC calculation
6) Bn*: the unnormalized virial (only for n>2)
7) error*: the error in the unnormalized virial

The data files contain:
1) index
2) density (in units for N/V)
3) alpha
4) S_nem
5) F(iso): free energy in the isotropic phase (in units of N*k_BT)
6) F(nem): free energy in the nematic phase (in units of N*k_BT)
7) P(iso): pressure in the isotropic phase (in units of N*k_BT)
8) P(nem): pressure in the nematic phase (in units of N*k_BT)
9) mu(iso): chemical potential in the isotropic phase (in units of N*k_BT)
10) mu(nem): chemical potential in the nematic phase (in units of N*k_BT)
