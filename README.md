# MMMC

The `python` codes in this repository are made available with https://arxiv.org/abs/XXXX.YYYYY

There are currently four different codes for matrix models discussed in the article above. 

1. `1MM.py` - One matrix model with quartic potential. If you want to change to cubic potential, you have to modify `def potential(X)`
and `def force(X)` appropriately. You can run this for 100 time units (also known as Molecular Dynamics (MD) time units)
as `python 1MM.py 0 1 300 100`

2. `2MM.py` - Hoppe-type models in the symmetric or symmetry broken phase. You can run an instance of
N = 300 model in symmetric phase for 10 time units by doing `python 2MM.py 0 1 300 10 1`

3. `MMC_3_4.py` - Matrix chain (open or closed) models with three or four matrices. If you want more than four matrices, you have to edit
the code. You can run an instance of N = 300 model (note that `NMAT = 3` is hardcoded, you can change it to `NMAT = 4`) 
for 10 time units with g = 1 and c = \kappa = 1.35 by doing `python MMC_3_4.py 0 1 300 10 1 1.35 1.35`

4. `YM_type.py` - This code can be used to study Yang-Mills type models with no mass terms. You can run an instance of
N = 300 model for 10 time units for D = 4 and \lambda = 1 by doing `python YM_type.py 0 1 300 10 4 1`


Some models like Yang-Mills with mass terms etc. are left for the interested reader. 
Please report any problem/bug in these codes by emailing at: raghav.govind.jha@gmail.com 


If you find this useful, please cite the following paper [[1]](#1). 


Last updated: 
Waterloo, 
October 21, 2021

## References
<a id="1">[1]</a> 
Jha, R. G. (2021). 


