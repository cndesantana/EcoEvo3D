Pkg.checkout("EcoEvo3D","model1")
using EcoEvo3D

function main()
	seed = 1;#seed for random numbers (to control the outputs)
	nreal = 10;#Number of realizations
	Gmax = 100;#Maximum number of Generations
	J = 10000;#Total Number of individuals in the system 
	anaG = round(Int64,J*2)::Int64;#Threshold to consider speciation (anagenesis and cladogenesis)
	v =  (0.00001 + rand()* 0.001)::Float64;#Cladogenesis speciation
	mr = (0.00001 + rand()* 0.0001)::Float64;#Regional migration - 
	ml = (0.002 + rand()* 0.2)::Float64;#Local Migration

	distmatfile = "upstreamcost_m3_0.txt";#Name of the file that contains the location of the points with clownfishes.
	verticesdata = "VerticesData.txt";#Name of the file that contains the size and the height of the points with clownfishes.
	model = 1;
		
	EcoEvo3D.dynamic(seed,nreal,Gmax,J,v,mr,ml,anaG,distmatfile,verticesdata,model);
end

main()
