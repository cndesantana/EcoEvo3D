using EcoEvo3D

function main()
	seed = 1;#seed for random numbers (to control the outputs)
	nreal = 1;#Number of realizations
	Gmax = 100;#Maximum number of Generations
	J = 10000;#Total Number of individuals in the system 
	distmatfile = "Dendritic_cost0_size100.txt";#Name of the file that contains the location of the points with clownfishes.
	verticesdata = "Dendritic_cost0_size100_elevation.txt";#Name of the file that contains the size and the height of the points with clownfishes.
	model = 1;
		
	EcoEvo3D.dynamic(seed,nreal,Gmax,J,distmatfile,verticesdata,model);
end

main()
