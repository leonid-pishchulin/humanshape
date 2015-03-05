HUMANSHAPE - 3D Human Shape Modeling, Fitting and Evaluation
=====

This short documentation describes steps necessary to compile and run the 3D human body shape building, manipulation, fitting and evaluation presented in the paper:

**Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobald and Bernt Schiele  
Building Statistical Shape Spaces for 3D Human Modeling
ArXiv, February 2015**

Compiling
---

This code was developed under Linux (Debian _wheezy_, 64 bit) and was tested only in this environment.  

1. Set `MATLAB_HOME` in following make files:

    ```
    external/lbfgsb-for-matlab/Makefile
    shapemodel/Makefile
    evaluation/statQuality/align.mk
    evaluation/statQuality/evaluation.mk
    ```
2. Switch to the top level directory of the source code, issue `make` from the command line

Getting the models
---
1. Download the models
```
    wget http://datasets.d2.mpi-inf.mpg.de/humanshape/caesar.zip
    wget http://datasets.d2.mpi-inf.mpg.de/humanshape/caesar-norm-wsx.zip
    wget http://datasets.d2.mpi-inf.mpg.de/humanshape/caesar-norm-nh.zip
```
2. Unzip the models
```
    unzip caesar.zip && rm -f caesar.zip
    unzip caesar-norm-wsx.zip && rm -f caesar-norm-wsx.zip
    unzip caesar-norm-nh.zip && rm -f caesar-norm-nh.zip
```

Running
---
1. Start matlab
2. Edit file `fitting/expParams.m`
   	1) point `p.rootDir` to the full path to the source code directory
   	2) point `p.modelInDir` to the model directory, e.g. `caesar/`
3. Run `demo`

TODO
---
1. Add evaluation code

   	1) per-vertex mean fitting accuracy

   	2) total fitting accuracy

   	3) compactness, generalization and specificity

If you have any questions, send an email to leonid@mpi-inf.mpg.de with a topic "humanshape".
