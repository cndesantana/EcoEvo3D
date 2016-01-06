#Pkg.checkout("EcoEvo3D","model1")
using EcoEvo3D

function main(ngenana_)
	seed = 177;#seed for random numbers (to control the outputs)
	nreal = 5;#Number of realizations
	Gmax = 10000;#Maximum number of Generations
	J = 10000;#Total Number of individuals in the system 
        nGenAna = ngenana_
	anaG = J*nGenAna;#Threshold to consider speciation (anagenesis and cladogenesis)
	v =  (0.000001 + rand()* 0.00001);#Cladogenesis speciation
	mr = (0.00001 + rand()* 0.0001);#Regional migration - 
	ml = (0.002 + rand()* 0.02);#Local Migration

	distmatfile = "upstream_cost0_size8_m3.txt";#Name of the file that contains the location of the points with clownfishes.
	verticesdata = "VerticesData.txt";#Name of the file that contains the size and the height of the points with clownfishes.
	model = 1;
		
	EcoEvo3D.dynamic(seed,nreal,Gmax,J,v,mr,ml,anaG,distmatfile,verticesdata,model);
end

if(length(ARGS)==0)
    println();
    println("Please give as an argment to the program the 'number of generations for anagenesis speciation':");
    println("To run: julia testm3.jl <anaG>");
    println();
else
    main(float(ARGS[1]))
end
