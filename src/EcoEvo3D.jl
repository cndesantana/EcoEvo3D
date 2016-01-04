# version 0.1

module EcoEvo3D 

function readLocations(filename)
	dat = readdlm(filename,'\t');
	@inbounds w = hcat(dat[:,2],dat[:,3]);
	w;
end

function readDistanceMatrix(filename)
	dat = readdlm(filename,' ');
	dat;
end

function changediagonal!(MA,nrows,value)
	@inbounds for(i in 1:nrows)
		@inbounds MA[i,i]=value;
	end
end

function createMRM(MRM,Sti,Sp,ts)
	MRM = [Sti Sp 1 ts];
	MRM;
end

function createMC(MC,Sti,Sp,ts)
	MC = [Sti Sp 1 ts];
	MC;
end

function createMA(MA,Sti,Stj,Sp,ts)
	MA = [Sti Stj Sp 1 ts];
	MA;
end

function updateMC(MC)
	if length(MC)>0
		MC[:,3] = MC[:,3].+1;#event: -1 for extinction; 1 for speciation
	end
	
	MC;
end

function updateMA(MA)
	if length(MA)>0
		MA[:,4]= MA[:,4] .+ 1;
	end
	MA;
end

function checkIfThereIsMRM(MRM,Sti,Sp,ts)
	if length(MRM)==0
		MRM = createMRM(MRM,Sti,Sp,ts);
	else
		pos = find( (MRM[:,1].==Sti) & (MRM[:,2].==Sp))
		if length(pos)==0
			MRM = cat(1,MRM,[Sti Sp 1 ts]);
		end
	end
	MRM;
end                                    

function checkIfThereIsMC(MC,Sti,Sp,ts)
	if length(MC)==0
		MC = createMC(MC,Sti,Sp,ts);
	else
		pos = find( (MC[:,1].==Sti) & (MC[:,2].==Sp))#position in the matrix MA referred to the presence of individuals of species 'Sp' coming from site 'Stj' to site 'Sti'
		if length(pos)==0
			MC = cat(1,MC,[Sti Sp 1 ts]);
		end
	end
	MC;
end

function printPhylogeny(new,old,ts,phylogenyfile,ri)
	writedlm(phylogenyfile,[ri old new ts],' '); 
end

function checkAna(MA,R,anaG,lastspecies,listofanagenesis,ts,phylogenyfile,ri)
	pos = [];
	Sti=0;
	if length(MA)>0
		pos = find(MA[:,4] .>= anaG)
		if length(pos)>0
		@inbounds for (a in 1:length(pos))
				@inbounds Sti=MA[pos[a],1];#pos represents the lines of the matrix. pos[a] is one line. MA[pos[a],1] is the a-th target site
				@inbounds Stj=MA[pos[a],2];#pos represents the lines of the matrix. pos[a] is one line. MA[pos[a],2] is the a-th source site 
				@inbounds Sp=MA[pos[a],3];#pos represents the lines of the matrix. pos[a] is one line. MA[pos[a],3] is the a-th species 
				MA, R, lastspecies,listofanagenesis = AnagenesisSpeciation(MA,R,Sti,Stj,Sp,lastspecies,listofanagenesis,ts,phylogenyfile,ri);#speciation in target site
			end
		end
	end
	Sti,MA,R,lastspecies,listofanagenesis,pos;
end

function UAB(MA,MC)
	MC = updateMC(MC); 
	MA = updateMA(MA);
	MC,MA;
end

function UAC(MA)
	MA = updateMA(MA);
	MA;	
end

function UARM(MA,MC)
	MC = updateMC(MC); 
	MA = updateMA(MA);
	MC,MA;
end

function UALM(MA,MC,R,Sti,Stj,Sp,anaG,retG,ts)
	MC = updateMC(MC); 
	MA = updateMA(MA);

	if length(MA)==0 #Si no hay MA
		MA = createMA(MA,Sti,Stj,Sp,ts);#Crea MA
	else#Si hay MA
		pos = find( (MA[:,1].==Sti) & (MA[:,2].==Stj) & (MA[:,3].==Sp))#position in the matrix MA referred to the presence of individuals of species 'Sp' coming from site 'Stj' to site 'Sti'
		if length(pos)==0 #No hay la linea
			indalive = length(find(R[Sti,:].==Sp))#Hay individuos de la specie sp vivos en el sitio Sti
			if (indalive == 0) #Checking if there are individuals of species 'Sp' alive in site 'Sti'
				MA = cat(1,MA,[Sti Stj Sp 1 ts]);#Crea la linea, empiezando con 1
			end
		else #Ya hay la linea
#			MA = MA[1:size(MA,1).!=pos,:];#Borra la linea 'pos' de la matriz MA!!
			MA[pos,4] = MA[pos,4] .- retG;#retards anagenesis by increasing the remaining 
#			MA = MA[:,1:size(MA,2).!=pos];#Borra la columna 'pos' de la matriz MA!!
		end
	end

	MC,MA,R;
end

########################

function AnagenesisSpeciation(MA,R,Sti,Stj,Sp,lastspecies,listofanagenesis,ts,phylogenyfile,ri)
	newspeciesAna = lastspecies + 1;#the id of the new species
	oldindividuals = find( R[Sti,:].==Sp )#the position of all the individuals of the 'old' species 'Sp' in the target site 
	R[Sti,oldindividuals] = newspeciesAna;#the speciation itself: all the individuals of former species 'Sp' in the target site are now from a new species 'newspeciesAna'
	printPhylogeny(newspeciesAna,Sp,ts,phylogenyfile,ri);
 
	pos = find( (MA[:,1].==Sti) & (MA[:,3].==Sp))#position in the matrix MA referred to the presence of individuals of species 'Sp' in site 'Sti' 
	MA = MA[1:size(MA,1).!=pos,:];#Borra la linea 'pos' de la matriz MA!!

	if length(listofanagenesis)>0
		listofanagenesis = cat(1,listofanagenesis,[Sti lastspecies]);
	else
		listofanagenesis = cat(1,[Sti lastspecies]);	
	end

	MA,R,newspeciesAna,listofanagenesis;
end

function CladogenesisEvent(MC,R,Sti,Individual,lastspecies,ts,phylogenyfile,ri)
	newspeciesClado = lastspecies + 1;
	printPhylogeny(newspeciesClado,R[Sti,Individual],ts,phylogenyfile,ri);
   	R[Sti,Individual] = newspeciesClado;
	MC = checkIfThereIsMC(MC,Sti,newspeciesClado,ts);
	MC,R,newspeciesClado; 
end

function LocalMigrationEvent(R,KillHab,MigrantHab,KillInd,Dc,Ji,S)
	MigrantSpecies = -1;

	allpos = find(Dc[KillHab,:] .>= MigrantHab);#All the sites at a distance lower than the threshold 'MigrantHab'
	kr = minimum(allpos[find(allpos .!= KillHab)]);
	MigrantInd = rand(1:Ji[kr]);
	MigrantSpecies = R[kr,MigrantInd];
	R[KillHab,KillInd] = R[kr,MigrantInd];
	return kr-1,R,MigrantSpecies;
end;

function RegionalMigrationEvent(MRM,R,Sti,Individual,ts,lastspecies)
	newspeciesMR = lastspecies + 1;
   	R[Sti,Individual] = newspeciesMR;
	MRM = checkIfThereIsMRM(MRM,Sti,newspeciesMR,ts);                                    

	return MRM,R,newspeciesMR; 
end;

function BirthEvent(R,BirthLocal,KillInd,KillHab)
	R[KillHab,KillInd] = R[KillHab,BirthLocal];
	R;
end

function calculateSpeciationMA(MA,listofanagenesis,R,S,Ji)
	speciatedMA = round(Int64,zeros(S));
	if (length(listofanagenesis) > 0)
		@inbounds for i in 1:S
			@inbounds speciesR = unique(sort(R'[1:Ji[i],i]))';
			pos = find(listofanagenesis[:,1].== i);
			if length(pos)>0
				speciesMA = unique(sort(listofanagenesis[pos,2]));
				speciatedMA[i] = length(setdiff(speciesMA,setdiff(speciesMA,speciesR))');
			end
		end
	end
#	writedlm("MAmatrix.dat",MA,' ');
	speciatedMA;
end


function calculateSpeciationMC(MC,R,S,k,Ji)
	speciatedMC = round(Int64,zeros(S));
	if (length(MC) > 0)
		@inbounds for i in 1:S
		@inbounds speciesR = unique(sort(R'[1:Ji[i],i]))';
			pos = find( (MC[:,1].== i) & (MC[:,3].>= k) );#Only the species that are in the system for more than k time steps (Protracted)
			if length(pos)>0
				speciesMC = unique(sort(MC[pos,2]));
				speciatedMC[i] = length(setdiff(speciesMC,setdiff(speciesMC,speciesR))');
			end
		end
	end
#	writedlm("MCmatrix.dat",MC,' ');
	speciatedMC;
end

function calculateSpeciationMR(MRM,R,S,Ji)
	speciatedMRM = round(Int64,zeros(S));
	if (length(MRM) > 0)
		@inbounds for i in 1:S
			@inbounds speciesR = unique(sort(R'[1:Ji[i],i]))';
			pos = find(MRM[:,1].== i);
			if length(pos)>0
				speciesMRM = unique(sort(MRM[pos,2]));
				@inbounds speciatedMRM[i] = length(setdiff(speciesMRM,setdiff(speciesMRM,speciesR))');
			end
		end
	end
#	writedlm("MRMmatrix.dat",MRM,' ');
	speciatedMRM;
end

function richnessanalysis!(S,R,Ji,richnessspeciesR,alpharich)	
	#%gamma richness
	@inbounds for (i in 1:S)
		@inbounds AR = sort(R'[1:Ji[i],i]);
		richnessspeciesR = [richnessspeciesR; unique(AR)];
		@inbounds alpharich[i] = length(unique(AR));
	end;
	gamma = length(unique(sort(richnessspeciesR)));
	gamma,alpharich;	
end

function dynamic(seed,nreal,Gmax,J,v,mr,ml,anaG,distmatfile,verticesdata,model)
	lastspecies = 0;
	DI = readDistanceMatrix(distmatfile);#the location of the points of the landscape.
	Dc = cumsum(DI,2);
	dT = sum(DI,2); 
	const S = length(DI[:,1]);#Number of sites
#### To get the cost
        substring = distmatfile[search(distmatfile,'_')+1:end];
        modelstop = search(substring,'_');
        coststop = search(substring,'t');
	const cost = float(substring[coststop+1:modelstop-1]);

	dat = readdlm(verticesdata,' ');#Proportion of individuals of each site
#	Pj = round(Int64,ones(S));#Proportion of individuals of each site - if an input file is not defined, the proportion is the same for each site
	Pj = dat[:,2];#In case we define different carrying capacities for each site
	mins = find(dat.==minimum(dat[:,1]));
	entrypoint = mins[rand(1:length(mins))];
	t=1;#Sites have different sizes and are located at different height.
	Ji=round(Int64,J * Pj/sum(Pj));	

	outputfile = open(string("RichnessPerSite_AnaG_",anaG,"_MR_",signif(mr,3),"_VR_",signif(v,3),".txt"),"a")	
	writedlm(outputfile,["Real Cost Model J G anaG retG Site Ji dT mr ml v gamma alpharich SpecANA SpecCLA SpecMR DispersalRich"]); 

	phylogenyfile = open(string("Phylogeny_AnaG_",anaG,"_MR_",signif(mr,3),"_VR_",signif(v,3),".txt"),"a")	
	writedlm(phylogenyfile,["Repl Ancestral Derived Age"]); 


	ts=0;

	R = round(Int64,zeros(S,maximum(Ji)));

	MA = Array(Int64,0,0);#Matrix to calculate Anagenesis speciation (MA = [Sti, Stj, S, C])
	MC = Array(Int64,0,0);#Matrix to control Cladogenesis speciation (MC = [Sti, S, E])
	MRM = Array(Int64,0,0);#Matrix to control events of Regional Migration (MRM = [Sti, S, ts])
	@inbounds for (ri in 1:nreal)#realizations
		srand(seed+(7*ri));
		s = rand();#prob to be a static landscape
		listofanagenesis = Array(Int64,0,0);
		
		#%Resources
#		G = rand(15:Gmax);#Minimum of 15 generations
		G = Gmax;#Minimum of 15 generations

		retG = round(Int64,anaG/1000);#gene flow retard in anagenetic speciation

	
		@inbounds for (k = 1:G)#%population-metapopulation-metacommunity dynamics (not-tracking multitrophic metacommunity dynamics!)
			ld = 0.0;# ld=0 because the landscape is static
			if ld > s;#%landscape dynamic: rgn
				#To change the cost 
			end#ld>s   

###### ANAGENESIS
			@inbounds for (j = 1:sum(Ji))#For each individual in each site
				ts = ts+1;#Time step increases
				mvb = rand();
		       		#Demography - Resources
		       		KillHab = rand(1:S);#which site to kill
				KillInd = rand(1:Ji[KillHab]);#which individual to kill
				MigrantHab = rand()*maximum(Dc[KillHab,:]);

				BirthLocal = rand(1:Ji[KillHab]);#which individual to born

		       		if mvb <= ml;#Local Migration event
					kr,R,MigrantSpecies = LocalMigrationEvent(R,KillHab,MigrantHab,KillInd,Dc,Ji,S);
					MC,MA,R = UALM(MA,MC,R,KillHab,kr,MigrantSpecies,anaG,retG,ts);#Update Anagenesis after Local Migration event
		       		elseif (mvb > ml) & (mvb <= ml+mr);#Regional Migration event
					MRM,R,lastspecies = RegionalMigrationEvent(MRM,R,entrypoint,KillInd,ts,lastspecies);#Speciation through the entry point
					MC,MA = UARM(MA,MC);#Update Anagenesis after Regional Migration event 
				elseif (mvb > ml+mr) & (mvb <= ml+mr+v)
					if(v > 0)#we only simulate Cladogenesis when the probability is higher than 0
						MC,R,lastspecies = CladogenesisEvent(MC,R,KillHab,KillInd,lastspecies,ts,phylogenyfile,ri);
						MA = UAC(MA);#Update Anagenesis after Cladogenesis Speciation
					end
		       		else #Birth event
					R = BirthEvent(R,BirthLocal,KillInd,KillHab);#Birth event
					MC,MA = UAB(MA,MC);#Update MA after Birth event
		       		end;
####ANAGENESIS
				lengthlistbefore = length(listofanagenesis);#Size of Anagenesis Matrix before checking for Anagenesis
				St,MA,R,lastspecies,listofanagenesis,pos = checkAna(MA,R,anaG,lastspecies,listofanagenesis,ts,phylogenyfile,ri); #After the update of matrix MA, we check for events of Anagenesis
				lengthlistafter = length(listofanagenesis);#Size of Anagenesis Matrix after checking for Anagenesis
####CLADOGENESIS
#				if ((lengthlistafter - lengthlistbefore) > 0)#If there is at least one Anagenesis Speciation Event
#					mvc = rand();
#					if (mvc <= v);#Probability to occur a Cladogenesis Speciation event
#						allInd = find(R[St,:].==lastspecies);#All individuals from the new species
#						if (length(allInd)>0)
#							KillInd = allInd[rand(1:length(allInd))];#Randomly choose an individual from the new species
#							MC,R,lastspecies = CladogenesisEvent(MC,R,KillHab,KillInd,lastspecies,ts,phylogenyfile,ri);#Cladogenesis event
#							MA = UAC(MA);#Update Anagenesis after Cladogenesis Speciation event
#						end #if allInd
#					end #if mvc
#				end #if lengthlistafter
			end;#end S*Ji
		end;#end Gmax
	
		#To analyze the resulting richness
		richnessspeciesR = Int64[];
		alpharich = zeros(S);
		gamma,alpharich = richnessanalysis!(S,R,Ji,richnessspeciesR,alpharich);
		SpecANA = calculateSpeciationMA(MA,listofanagenesis,R,S,Ji);
		SpecCLA = calculateSpeciationMC(MC,R,S,anaG,Ji);
		SpecMR = calculateSpeciationMR(MRM,R,S,Ji);
		DispersalRich = alpharich - (SpecANA + SpecCLA + SpecMR);
		for i in 1:S
			writedlm(outputfile,[ri cost model J G anaG retG i Ji[i] dT[i] mr ml v gamma alpharich[i] SpecANA[i] SpecCLA[i] SpecMR[i] DispersalRich[i]],' '); 
			flush(outputfile);
		end 
	end#%ri
	close(outputfile);
#	close(logfile);
end#ecoevo3d function

end #module
