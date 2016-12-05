This is a ** two way Fluid Structure interaction (FSI) coupling code** based on *immersed boudnary method*. It it a good start to know FSI.
Its main purpose is to study blood flow in microcirculation.
The fluid field is solved through *Lattice Boltzmann method*, while the cells 
are modeled as a *membrane network consisting of beads and springs*. Two unique
motions related to cells under shear, namely tumbling and tank treading,
can be observed. The effect of RBCs on particle dispersion rate were studied 
and published. If you feel this code is useful for your research, please consider
cite this paper:

Tan, Jifu, et al. [Characterization of nanoparticle dispersion in red blood cell suspension by the lattice boltzmann-immersed boundary method!](http://www.mdpi.com/2079-4991/6/2/30/htm) Nanomaterials 6.2 (2016): 30.

#Requirements:
Linux OS, g++ compiler, matlab for pre/postprocessing inputs.

It is preferred running the code in linux computers. I haven't tested it in windows. A little extrawork may be required to run it in windows. I believe it should be straightforward. 
To fully understand the code, you should know c++, familiar with Lattice Boltzmann method, immersed boundary method, spring-bead models(stretching, bending) such that you can check the source code and know what it does. Then you can modify the code whatever you need. The above reference should give you a good introduction and background to what it does in the code.  

#How to run it?
Download the source file and type `make` in the terminal, it should compile. Next, just type `./fsi` in the terminal. It should run.  

It requires input files: Velchannel.txt, MultiCells.txt particles.txt
Those are input files for fluid channel, 2D red blood cells, and point particels. These files can be created using the matlab scripts in the input/postprocessing folder. What it does is to simulate the particle cell mixture in a straight channel, as shown in the above mentioned paper. You may need to change it based on what you need. 

You also need to change the *main.cpp* file based on how many solid components you want to include in your simulation. Here we show two components, cells and particles. You can include red blood cells, White blood cells, for example. 

Notice it is a research code for my study only, thus, it may not provide enough documentation or tutorials. However, I would like to help if you find any problems/bugs, such as segmentation fault. Shoot me email if you have any questions, typically I would get back to you in 24 hours.  

