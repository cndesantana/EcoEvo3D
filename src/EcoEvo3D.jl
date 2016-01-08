# version 0.1

module EcoEvo3D

using StatsBase

function getSampleFromArray(items)
  return sample(items)
end

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
	return MC;
end

function createMA(MA,Sti,Stj,Sp,anaG)
	MA = [Sti Stj Sp anaG];
	return MA;
end

function updateMC(MC)
	if length(MC)>0
		MC[:,3] = MC[:,3].+1;#event: -1 for extinction; 1 for speciation
	end
	return MC;
end

function updateMA(MA)
	if length(MA)>0
		MA[:,4] = MA[:,4] .- 1;
	end
	return MA;
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
	return MRM;
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
	return MC;
end

function GetGammaRichness(R,S)
	richnessspeciesR = [];
	#%gamma richness
	for (i in 1:S)
		AR = sort(R[i])
		richnessspeciesR = [richnessspeciesR; unique(AR)];
	end;
	extantR = length(unique(sort(richnessspeciesR)));
	gamma = extantR;

	return gamma;
end

function OutputPerGeneration(outputfilepergen,ri,cost,J,G,S,k,anaG,retG,mr,ml,v,gamma,alpharich,SpecANA,SpecCLA,SpecMR,DispersalRich)
	for i in 1:S
		writedlm(outputfilepergen, [ri cost J G k anaG retG mr ml v gamma i alpharich[i] SpecANA[i] SpecCLA[i] SpecMR[i] DispersalRich[i]],' ');
	end
	flush(outputfilepergen);#To print in the output file each realization
	return;
end

function printPhylogeny(phylogenyfile,old,new,ts,ri)
	writedlm(phylogenyfile,[ri old new ts],' ');
        flush(phylogenyfile);
end

function UAB(MA,MC)
	MC = updateMC(MC);
	MA = updateMA(MA);
	return MC,MA;
end

function UAC(MA)
	MA = updateMA(MA);
	return MA;
end

function UARM(MA,MC)
	MC = updateMC(MC);
	MA = updateMA(MA);
	return MC,MA;
end

function UALM(MA,MC,R,Sti,Stj,Sp,anaG,retG,ts)
	MC = updateMC(MC);
	MA = updateMA(MA);

	if length(MA)==0 #Si no hay MA
  	MA = createMA(MA,Sti,Stj,Sp,anaG);#Create MA
	else#Si hay MA
		pos = find( (MA[:,1].==Sti) & (MA[:,2].==Stj) & (MA[:,3].==Sp))#position in the matrix MA referred to the presence of individuals of species 'Sp' coming from site 'Stj' to site 'Sti'
		if length(pos)==0 #No hay la linea
			indalive = length(find(R[Sti].==Sp))#Hay individuos de la specie sp vivos en el sitio Sti
			if (indalive == 0) #Checking if there are individuals of species 'Sp' alive in site 'Sti'
#        			println("<<<<<<<<<<< A D D    R O W   T O    A N A G E N E S I S   <<<<<<<<<<<<<<<<<");
				MA = cat(1,MA,[Sti Stj Sp anaG]);#Create the row in MA matrix, starting by 1
      			end
		else #Ya hay la linea
#			println("<<<<<<<<<<< R E T A R D    A N A G E N E S I S   <<<<<<<<<<<<<<<<<");
			MA[pos,4] = MA[pos,4] .+ retG;#retards anagenesis by increasing the remaining
		end
	end
  	return MC,MA,R;
end

########################

function checkAna(MA,R,anaG,lastspecies,listofanagenesis,ts,phylogenyfile,ri)
	pos = [];
	if length(MA)>0
		pos = find(MA[:,4] .<= 0)#is there any
		if length(pos)>0
		  @inbounds for (p in pos)
  			Sti=round(Integer,MA[p,1]);#pos represents the rows of the matrix. p is one row. MA[p,1] is the target site
  			Stj=round(Integer,MA[p,2]);#pos represents the rows of the matrix. p is one row. MA[p,2] is the source site
		  	Sp=round(Integer,MA[p,3]);#pos represents the rows of the matrix. p is one row. MA[p,3] is the ancient species
			  MA, R, lastspecies,listofanagenesis = AnagenesisSpeciation(MA,R,Sti,Stj,Sp,lastspecies,listofanagenesis,ts,phylogenyfile,ri);#speciation in target site
#			  MA, R, lastspecies,listofanagenesis = MA, R, lastspecies, listofanagenesis;
		  end
	  end
  end
  return MA,R,lastspecies,listofanagenesis;
end

function AnagenesisSpeciation(MA,R,Sti,Stj,Sp,lastspecies,listofanagenesis,ts,phylogenyfile,ri)
	newspeciesAna = lastspecies + 1;#the id of the new species
	ancientindividuals = find( R[Sti].==Sp )#the position of all the individuals of the ancient species 'Sp' in the target site
	R[Sti][ancientindividuals] = newspeciesAna;#the speciation itself: all the individuals of former species 'Sp' in the target site are now from a new species 'newspeciesAna'
        println("new ANAgenesis Speciation : Species ",newspeciesAna," - time: ",ts);
 	pos = find( (MA[:,1].==Sti) & (MA[:,3].==Sp))#position in the matrix MA referred to the presence of individuals of species 'Sp' in site 'Sti'
	MA = MA[1:size(MA,1).!=pos,:];#Borra la linea 'pos' de la matriz MA!!
	printPhylogeny(phylogenyfile,Sp,newspeciesAna,ts,ri);

  if length(listofanagenesis)>0
		listofanagenesis = cat(1,listofanagenesis,[Sti lastspecies]);
	else
		listofanagenesis = cat(1,[Sti lastspecies]);
	end

	return MA,R,newspeciesAna,listofanagenesis;
end

function CladogenesisEvent(MC,R,Sti,Individual,lastspecies,ts,phylogenyfile,ri)
	newspeciesClado = lastspecies + 1;
	printPhylogeny(phylogenyfile,R[Sti][Individual],newspeciesClado,ts,ri);
        println("new CLADOgenesis Speciation : Species ",newspeciesClado," - time: ",ts);
   	R[Sti][Individual] = newspeciesClado;
	MC = checkIfThereIsMC(MC,Sti,newspeciesClado,ts);
	return MC,R,newspeciesClado;
end

function LocalMigrationEvent(R,KillHab,MigrantHab,KillInd,Dc,Ji,S)
	MigrantSpecies = -1;

	allpos = find(Dc[KillHab,:] .>= MigrantHab);#All the sites at a distance lower than the threshold 'MigrantHab'
	kr = minimum(allpos[find(allpos .!= KillHab)]);
	MigrantInd = rand(1:Ji[kr]);
	MigrantSpecies = R[kr][MigrantInd];
#        println("LOCAL Migration");
	R[KillHab][KillInd] = R[kr][MigrantInd];
	return kr,R,MigrantSpecies;
end;

function RegionalMigrationEvent(MRM,R,Sti,Individual,ts,lastspecies)
  newspeciesMR = lastspecies + 1;
  R[Sti][Individual] = newspeciesMR;
  println("new RegionalMig Speciation: Species ",newspeciesMR," - time: ",ts)
  MRM = checkIfThereIsMRM(MRM,Sti,newspeciesMR,ts);

  return MRM,R,newspeciesMR;
end;

function BirthEvent(R,BirthLocal,KillInd,KillHab)
	R[KillHab][KillInd] = R[KillHab][BirthLocal];
	return R;
end

function calculateSpeciationMA(MA,listofanagenesis,R,S,Ji)
	speciatedMA = zeros(S);
	if (length(listofanagenesis) > 0)
		@inbounds for i in 1:S
			@inbounds speciesR = unique(sort(R[i]))';
			pos = find(listofanagenesis[:,1].== i);
			if length(pos)>0
				speciesMA = unique(sort(listofanagenesis[pos,2]));
				speciatedMA[i] = length(setdiff(speciesMA,setdiff(speciesMA,speciesR))');
			end
		end
	end
	return speciatedMA;
end


function calculateSpeciationMC(MC,R,S,k,Ji)
	speciatedMC = zeros(S);
	if (length(MC) > 0)
		@inbounds for i in 1:S
		@inbounds speciesR = unique(sort(R[i]))';
			pos = find( (MC[:,1].== i) & (MC[:,3].>= k) );#Only the species that are in the system for more than k time steps (Protracted)
			if length(pos)>0
				speciesMC = unique(sort(MC[pos,2]));
				speciatedMC[i] = length(setdiff(speciesMC,setdiff(speciesMC,speciesR))');
			end
		end
	end
	return speciatedMC;
end

function calculateSpeciationMR(MRM,R,S,Ji)
	speciatedMRM = zeros(S);
	if (length(MRM) > 0)
		@inbounds for i in 1:S
      @inbounds speciesR = unique(sort(R[i]))';
      pos = find(MRM[:,1].== i);#check for site i
      if length(pos)>0
        speciesMRM = unique(sort(MRM[pos,2]));
        @inbounds speciatedMRM[i] = length(setdiff(speciesMRM,setdiff(speciesMRM,speciesR))');
      end
		end
	end
#	writedlm("MRMmatrix.dat",MRM,' ');
	return speciatedMRM;
end

function richnessanalysis!(S,R,Ji,richnessspeciesR,alpharich)
	#%gamma richness
	@inbounds for (i in 1:S)
		@inbounds AR = sort(R[i]);
		richnessspeciesR = [richnessspeciesR; unique(AR)];
		@inbounds alpharich[i] = length(unique(AR));
	end;
	gamma = length(unique(sort(richnessspeciesR)));
	return gamma,alpharich;
end

function normalizeDc!(Dc,S)
  for(i in 1:S)
    Dc[i,:] = Dc[i,:]/maximum(Dc[i,:]);
  end
end

function repeatval(x,N)
  return map(v -> v=x,1:N);
end

#Each site starts with one different species
function initialPopulation(R,Ji)
  lastspecies = 0;
  for(iJ in Ji)#iterating over columns
     lastspecies = lastspecies+1;
     popi = repeatval(lastspecies,iJ);
     push!(R,popi);
  end
  return R,lastspecies;
end

function dynamic(seed,nreal,Gmax,J,v,mr,ml,anaG,retG,distmatfile,verticesdata,model)
	DI = readDistanceMatrix(distmatfile);#the location of the points of the landscape.
	Dc = cumsum(DI,2);#Should we divide
  dT = sum(DI,2);
	const S = length(DI[:,1]);#Number of sites
  normalizeDc!(Dc,S);
#### To get the cost
        substring = distmatfile[search(distmatfile,'_')+1:end];
        modelstop = search(substring,'_');
        coststop = search(substring,'t');
	const cost = float(substring[coststop+1:modelstop-1]);

	sitesdata = readdlm(verticesdata,' ',header=true)[1];#Proportion of individuals of each site
#	Pj = ones(S);#Proportion of individuals of each site - if an input file is not defined, the proportion is the same for each site
	Pj = sitesdata[:,2];#In case we define different carrying capacities for each site we consider the volume of the lakes (third column)
	mins = find(sitesdata.==minimum(sitesdata[:,1]));#the first column represents the height of the lakes
	entrypoint = mins[rand(1:length(mins))];
	t=1;#Sites have different sizes and are located at different height.
	Ji=round(Integer,J * Pj/sum(Pj));

	outputfilepergen = open(string("RichnessPerGen_AnaG_",anaG,"_cost_",cost,"_MR_",signif(mr,3),"_VR_",signif(v,3),".txt"),"a")
	writedlm(outputfilepergen, ["Real Cost J G Gi anaG retG mr ml v gamma Site alpharich SpecANA SpecCLA SpecMR DispersalRich"]);

	outputfile = open(string("RichnessPerSite_AnaG_",anaG,"_cost_",cost,"_MR_",signif(mr,3),"_VR_",signif(v,3),".txt"),"a")
	writedlm(outputfile,["Real Cost Model J G anaG retG Site Ji dT mr ml v gamma alpharich SpecANA SpecCLA SpecMR DispersalRich"]);

	phylogenyfile = open(string("Phylogeny_AnaG_",anaG,"_cost_",cost,"_MR_",signif(mr,3),"_VR_",signif(v,3),".txt"),"a")
	writedlm(phylogenyfile,["Repl Ancestral Derived Age"]);

  @inbounds for (ri in 1:nreal)#realizations
 	  lastspecies = 0;
  	ts=0;
    R=[];
    R,lastspecies = initialPopulation(R,Ji)
	  MA = Array(Number,0,0);#Matrix to calculate Anagenesis speciation (MA = [Sti, Stj, S, C])
	  MC = Array(Number,0,0);#Matrix to control Cladogenesis speciation (MC = [Sti, S, E])
	  MRM = Array(Number,0,0);#Matrix to control events of Regional Migration (MRM = [Sti, S, ts])
		srand(seed+(7*ri));
		listofanagenesis = Array(Number,0,0);

		#%Resources
#		G = rand(15:Gmax);#Minimum of 15 generations
		G = Gmax;#Minimum of 15 generations

		@inbounds for (k = 1:G)#%population-metapopulation-metacommunity dynamics (not-tracking multitrophic metacommunity dynamics!)
			@inbounds for (j = 1:sum(Ji))#For each individual in each site
				ts = ts+1;#Time step increases
				mvb = rand();
		       		#Demography - Resources
    		KillHab = rand(1:S);#which site to kill
				KillInd = rand(1:Ji[KillHab]);#which individual to kill
        KillIndEntryPoint = rand(1:Ji[entrypoint]);
				MigrantHab = rand()*maximum(Dc[KillHab,:]);

				BirthLocal = getSampleFromArray(1:length(R[KillHab]));#which individual to born

		    if mvb <= ml;#Local Migration event
					kr,R,MigrantSpecies = LocalMigrationEvent(R,KillHab,MigrantHab,KillInd,Dc,Ji,S);
					MC,MA,R = UALM(MA,MC,R,KillHab,kr,MigrantSpecies,anaG,retG,ts);#Update Anagenesis after Local Migration event
		    elseif (ml < mvb <= (ml+mr));#Regional Migration event
					MRM,R,lastspecies = RegionalMigrationEvent(MRM,R,entrypoint,KillIndEntryPoint,ts,lastspecies);#Speciation through the entry point
					MC,MA = UARM(MA,MC);#Update Anagenesis after Regional Migration event
				elseif ((ml+mr) < mvb <= (ml+mr+v))
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
				MA,R,lastspecies,listofanagenesis = checkAna(MA,R,anaG,lastspecies,listofanagenesis,ts,phylogenyfile,ri); #After the update of matrix MA, we check for events of Anagenesis
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
			richnessspeciesR = [];
			alpharich = zeros(S);
			gamma,alpharich = richnessanalysis!(S,R,Ji,richnessspeciesR,alpharich);
			SpecANA = calculateSpeciationMA(MA,listofanagenesis,R,S,Ji);
			SpecCLA = calculateSpeciationMC(MC,R,S,anaG,Ji);
			SpecMR = calculateSpeciationMR(MRM,R,S,Ji);
			DispersalRich = alpharich - (SpecANA + SpecCLA + SpecMR);
			OutputPerGeneration(outputfilepergen,ri,cost,J,G,S,k,anaG,retG,mr,ml,v,gamma,alpharich,SpecANA,SpecCLA,SpecMR,DispersalRich);
		end;#end Gmax

		#To analyze the resulting richness
		richnessspeciesR = [];
		alpharich = zeros(S);
		gamma,alpharich = richnessanalysis!(S,R,Ji,richnessspeciesR,alpharich);
		SpecANA = calculateSpeciationMA(MA,listofanagenesis,R,S,Ji);
		SpecCLA = calculateSpeciationMC(MC,R,S,anaG,Ji);
		SpecMR = calculateSpeciationMR(MRM,R,S,Ji);
		DispersalRich = alpharich - (SpecANA + SpecCLA + SpecMR);
		for i in 1:S
			writedlm(outputfile,[ri cost model J G anaG retG i Ji[i] dT[i] mr ml v gamma alpharich[i] SpecANA[i] SpecCLA[i] SpecMR[i] DispersalRich[i]],' ');
		end
		flush(outputfile);
	end#%ri
	close(outputfile);
	close(outputfilepergen);
	close(phylogenyfile);
#	close(logfile);
end#ecoevo3d function

end #module
