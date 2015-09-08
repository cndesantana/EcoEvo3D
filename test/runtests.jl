using EcoEvo3D

function main()
	seed = 1;#seed for random numbers (to control the outputs)
	nreal = 1;#Number of realizations
	Gmax = 100;#Maximum number of Generations
	J = 10000;#Total Number of individuals in the system 
	anaG = round(Int64,J*2)::Int64;#Threshold to consider speciation (anagenesis and cladogenesis)
	v =  (0.0001 + rand()* 0.001)::Float64;
	mr = (0.00001 + rand()* 0.0001)::Float64 
	ml = (0.02 + rand()* 0.2)::Float64;

	distmatfile = "Dendritic_cost0_size100.txt";#Name of the file that contains the location of the points with clownfishes.
	verticesdata = "Dendritic_cost0_size100_elevation.txt";#Name of the file that contains the size and the height of the points with clownfishes.
	model = 1;
		
	EcoEvo3D.dynamic(seed,nreal,Gmax,J,v,mr,ml,anaG,distmatfile,verticesdata,model);
end

main()
