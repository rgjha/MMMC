# MMMC

The `python` codes in this repository are made available with https://arxiv.org/abs/XXXX.YYYYY

There are currently four different codes for matrix models discussed in the article above in detail. 
The general form of the partition function is: <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{80}&space;\bg_black&space;\fn_phv&space;Z&space;=&space;\int&space;dM&space;\exp\Big[-N&space;\mbox{Tr}&space;V(M)\Big]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{80}&space;\bg_black&space;\fn_phv&space;Z&space;=&space;\int&space;dM&space;\exp\Big[-N&space;\mbox{Tr}&space;V(M)\Big]" title="Z = \int dM \exp\Big[-N \mbox{Tr} V(M)\Big]" /></a>
where `V` is the potential and `N` is the size of the Hermitian matrix. 


1. `1MM.py` - One matrix model with quartic potential defined as: <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{80}&space;\bg_black&space;\fn_phv&space;\large&space;V(M)&space;=&space;M^2/2&space;&plus;&space;gM^4/4" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{80}&space;\bg_black&space;\fn_phv&space;\large&space;V(M)&space;=&space;M^2/2&space;&plus;&space;gM^4/4" title="\large V(M) = M^2/2 + gM^4/4" /></a>. 
If you want to change to cubic potential, you have to modify `def potential(X)`
and `def force(X)` appropriately. You can run this for 100 time units (also known as Molecular Dynamics (MD) time units)
as `python 1MM.py 0 1 300 100`

2. `2MM.py` - Hoppe-type models in the symmetric or symmetry broken phase give by the potetial: <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{80}&space;\bg_black&space;\fn_phv&space;\large&space;V(X,Y)&space;=&space;X^2&space;&plus;&space;Y^2&space;-&space;h^2&space;[X,Y]^2&space;&plus;&space;gX^4&space;&plus;&space;gY^4" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{80}&space;\bg_black&space;\fn_phv&space;\large&space;V(X,Y)&space;=&space;X^2&space;&plus;&space;Y^2&space;-&space;h^2&space;[X,Y]^2&space;&plus;&space;gX^4&space;&plus;&space;gY^4" title="\large V(X,Y) = X^2 + Y^2 - h^2 [X,Y]^2 + gX^4 + gY^4" /></a>. 
You can run an instance of `N = 300` model in symmetric phase for 10 time units by doing `python 2MM.py 0 1 300 10 1`

3. `MMC_3_4.py` - Matrix chain (open or closed) models with three or four matrices. The potential is given by (for `p=4`):
<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{80}&space;\bg_black&space;\fn_phv&space;\large&space;V(M_{1},M_{2},M_{3},M_{4})&space;=&space;\Bigg(\sum_{i=1}^{4}&space;-M_{i}^2&space;-&space;g&space;M_{i}^{4}&space;&plus;&space;c&space;\sum_{i=1}^{3}&space;M_{i}M_{i&plus;1}&space;&plus;&space;\kappa&space;M_{4}M_{1}&space;\Bigg)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{80}&space;\bg_black&space;\fn_phv&space;\large&space;V(M_{1},M_{2},M_{3},M_{4})&space;=&space;\Bigg(\sum_{i=1}^{4}&space;-M_{i}^2&space;-&space;g&space;M_{i}^{4}&space;&plus;&space;c&space;\sum_{i=1}^{3}&space;M_{i}M_{i&plus;1}&space;&plus;&space;\kappa&space;M_{4}M_{1}&space;\Bigg)" title="\large V(M_{1},M_{2},M_{3},M_{4}) = \Bigg(\sum_{i=1}^{4} -M_{i}^2 - g M_{i}^{4} + c \sum_{i=1}^{3} M_{i}M_{i+1} + \kappa M_{4}M_{1} \Bigg)" /></a>. If you want more than four matrices, you have to edit the code. You can run an instance of N = 300 model (note that `NMAT = 3` is hardcoded, you can change it to `NMAT = 4`) 
for 10 time units with `g = 1` and `c =  κ = 1.35` by doing `python MMC_3_4.py 0 1 300 10 1 1.35 1.35`

4. `YM_type.py` - This code can be used to study Yang-Mills type models (with `D` matrices) with potential:
<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{80}&space;\bg_black&space;\fn_phv&space;\large&space;V&space;=&space;\Bigg(&space;\sum_{i&space;\neq&space;j}[X_i,&space;X_j]^{2}\Bigg)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{80}&space;\bg_black&space;\fn_phv&space;\large&space;V&space;=&space;\Bigg(&space;\sum_{i&space;\neq&space;j}[X_i,&space;X_j]^{2}\Bigg)" title="\large V = \Bigg( \sum_{i \neq j}[X_i, X_j]^{2}\Bigg)" /></a>. You can run an instance of `N = 300` model for 10 time units for` D = 4` and `λ  = 1` by doing `python YM_type.py 0 1 300 10 4 1`


Some models like Yang-Mills with mass terms etc. are left for the interested reader. 
Please write to me about any bug or comments: raghav.govind.jha@gmail.com 


If you find this repository useful, please cite the following paper [[1]](#1). 


Last updated: 
Waterloo, 
October 21, 2021


## References
<a id="1">[1]</a> 
Jha, R. G. (2021). 
