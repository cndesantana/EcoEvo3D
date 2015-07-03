# EcoEvo3D

To install EcoEvo3D in your machine:

1) Open a Julia terminal

2) run the following command: Pkg.clone("https://github.com/cndesantana/EcoEvo3D.git")

It will install the model in your default Julia folder (typically it will be in the ~/.julia/v0.X/EcoEvo3D folder - where 0.X is your current version of Julia)

3) You will find a test script in the 'tests' folder called runtests.jl. You can copy and paste the commands that are in the runtests.jl file in your julia terminal.

4) After defining the values of the parameters, just call the function EcoEvo3D.dynamic giving the parameters you want. 

5) Please note that in the 'tests' folder there are also 2 input files, Dendritic_cost0_size100_elevation.txt and Dendritic_cost0_size100.txt. You can use them to test the model.

