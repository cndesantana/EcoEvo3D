#Pkg.checkout("EcoEvo3D","model1")
using EcoEvo3D

function main(ngenana_,ngenretana_)
	seed = 177;#seed for random numbers (to control the outputs)
	nreal = 2;#Number of realizations
	Gmax = 2000;#Maximum number of Generations
	J = 50000;#Total Number of individuals in the system 
        nGenAna = ngenana_
	anaG = J*nGenAna;#Threshold to consider speciation (anagenesis and cladogenesis)
        nGenRetG = ngenretana_; #parameter to define gene flow retard in anagenetic speciation (10% of retard, minimum value of retard is 1)
	retG = nGenRetG*J;#gene flow retard in anagenetic speciation

	srand(seed);
	v =  (0.000005 + rand()* 0.00005);#Cladogenesis speciation
	mr = (0.00001 + rand()* 0.0001);#Regional migration - 
	ml = (0.008 + rand()* 0.008);#Local Migration

	distmatfile = "upstream_cost0.001_size8_m1.txt";#Name of the file that contains the location of the points with clownfishes.
	verticesdata = "VerticesData.txt";#Name of the file that contains the size and the height of the points with clownfishes.
#	verticesdata = "vert.txt";#Name of the file that contains the size and the height of the points with clownfishes.
	model = 1;
		
	EcoEvo3D.dynamic(seed,nreal,Gmax,J,v,mr,ml,anaG,retG,distmatfile,verticesdata,model);
end

if(length(ARGS)==0)
    println();
    println("Please give the correct parameters:");
    println("To run: julia test.jl <ngenana> <ngenretana>");
    println();
    println("where: ");
    println("- ngenana: number of generations for anagenesis speciation");
    println("- ngenretana: number of generations in the retard of an anagenesis speciation");
    println();
else
    main(float(ARGS[1]),float(ARGS[2]))
end
