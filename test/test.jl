#Pkg.checkout("EcoEvo3D","model1")
using EcoEvo3D

function main(ngenana_)
	seed = 177;#seed for random numbers (to control the outputs)
	nreal = 1;#Number of realizations
	Gmax = 100000;#Maximum number of Generations
	J = 10000;#Total Number of individuals in the system 
        nGenAna = ngenana_
	anaG = J*nGenAna;#Threshold to consider speciation (anagenesis and cladogenesis)
	v =  (0.0005 + rand()* 0.00005);#Cladogenesis speciation
	mr = (0.001 + rand()* 0.0001);#Regional migration - 
	ml = (0.08 + rand()* 0.008);#Local Migration

	distmatfile = "upstream_cost0.001_size8_m1.txt";#Name of the file that contains the location of the points with clownfishes.
	verticesdata = "VerticesData.txt";#Name of the file that contains the size and the height of the points with clownfishes.
#	verticesdata = "vert.txt";#Name of the file that contains the size and the height of the points with clownfishes.
	model = 1;
		
	EcoEvo3D.dynamic(seed,nreal,Gmax,J,v,mr,ml,anaG,distmatfile,verticesdata,model);
end

if(length(ARGS)==0)
    println();
    println("Please give as an argment to the program the 'number of generations for anagenesis speciation':");
    println("To run: julia test.jl <anaG>");
    println();
else
    main(float(ARGS[1]))
end
