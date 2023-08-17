# MMMC

There are four codes in this repository that accompany the article `Introduction to Monte Carlo for Matrix Models` 
available at https://arxiv.org/abs/2111.02410 or at https://scipost.org/10.21468/SciPostPhysLectNotes.46
 
The general form of the partition function is: 
 ```math
Z = \int dM \exp\Big[-N \text{Tr} V(M)\Big],
```
where `V` is the potential and `N` is the size of the Hermitian matrix. 

To use the code, you need Python with the required libraries. If you don't have these already, then it is easiest 
to just install Anaconda which is a cross-platform Python distribution for scientific computing
and available here: https://docs.anaconda.com/anaconda/install/
It usually works well and should be sufficient for all the codes in this repository. 

1. `1MM.py` - One-matrix (Hermitian) model with quartic potential defined as:

 ```math
V(M) = \frac{1}{2}M^2 + \frac{g}{4} M^4,
```
If you want to change to cubic potential, you have to modify `def potential(X)`
and `def force(X)` appropriately. You can run this for 100 time units (also known as Molecular Dynamics (MD) time units) with  N = 300
as `python 1MM.py 0 1 300 100`. See the article for more details. 

2. `2MM.py` - Hoppe-type models in the symmetric or symmetry broken phase give by the potential: 

 ```math
V(X,Y) = X^2 + Y^2 - h^{2}[X,Y]^{2} + gX^4 + gY^4,
```

You can run an instance of N = 150 model in symmetric phase for 10 time units by doing `python 2MM.py 0 1 150 10 1`. The last argument specifies whether we want to consider symmetric phase (1) or not (0). 

3. `3_4MMC.py` - Matrix chain (open or closed) models with three or four matrices. The potential is given by (for `p=4`):

 ```math
V(M_{1},M_{2},M_{3},M_{4}) = \Bigg(\sum_{i=1}^{4} -M_{i}^2 - g M_{i}^{4} + c \sum_{i=1}^{3} M_{i}M_{i+1} + \kappa M_{4}M_{1} \Bigg)
```


If you want more than four matrices, you have to edit the code. You can run an instance of N = 100 model (note that `NMAT = 3` is hard-coded, you can change it to `NMAT = 4`) for 10 time units with `g = 1` and `c =  κ = 1.35` by doing `python 3_4MMC.py 0 1 100 10 1 1.35 1.35`

4. `YMtype.py` - This code can be used to study Yang-Mills type models (with `D` matrices) with commutator potential. You can run an instance of N = 100 model for 10 time units for` D = 4` and `λ  = 1` by doing `python YMtype.py 0 1 100 10 4 1`


Some other models like Yang-Mills with mass terms etc. are left for the interested reader. 


## Contact, Cite et al. 

Please write to raghav.govind.jha@gmail.com for bug reports and comments. If you find this repository useful, please consider citing the paper [[1]](#1) and give this repository a star! We paste the bibtex entry below for you to copy if needed. 

```bibtex
 @article{Jha:2021exo,
    author = "Jha, Raghav G.",
    title = "{Introduction to Monte Carlo for Matrix Models}",
    eprint = "2111.02410",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    doi = "10.21468/SciPostPhysLectNotes.46",
    month = "11",
    year = "2021"
}
```


## About this work

This work started when the author was invited to give set of lectures on some numerical topic at a summer conference organized at Rensselaer Polytechnic institute (RPI). The author decided to talk on "Numerical aspects of matrix models". After few months, the author noticed a new paper by Volodya Kazakov & Zechuan Zheng: https://arxiv.org/abs/2108.04830. They applied numerical bootstrap to an unsolved matrix model and claimed several digits of precision for the moment of matrices. Volodya gave a talk at Perimeter and after the talk, I pulled my lecture codes and studied this model and found very good agreement with their results. I showed this to Pedro Vieira and he suggested that it would be nice to check other results in that paper and then write an introduction giving all the codes which I used so that future bootstrappers could check their result readily if required. You can see the reference of this work in the talk starting at 1:05:36 by Kazakov, available here: https://www.youtube.com/watch?v=34z3xE11ycY&t=3936s&ab_channel=BootstrapCollaboration about two weeks before the paper was uploaded on arXiv. 



## References
<a id="1">[1]</a> 
R. G. Jha, Introduction to Monte Carlo for Matrix Models, SciPost Physics Lecture Notes 46 (2022), arXiv:2111.02410 



Last updated: 
November 04, 2021
