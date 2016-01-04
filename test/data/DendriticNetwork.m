%------------------------------------------------------------------------
%Melian, C. J. @EAWAG, 2013
%Extension of Island Biogeography theory to dendritic networks
%Extension Rosindel & Phillimore (2012)
%1st: March-April 2013 2nd Jan-Feb 2014
%------------------------------------------------------------------------

%1. Input data: pairwise distances, number of species in each site...
HydroDistanceMatrix;%generate Migration-Hydrological-Matrix M

%2.Initial conditions
rand('seed',sum(100*clock));
nruns = 1000;%Number of replicates
SP = size(M);S = SP(1,1);%number of lakes
J = 10000;%number of individuals per lake
TG = 100000;%generations

%3. Metacommunity dynamics
for rr = 1:nruns;
    G = unidrnd(TG,1,1);%Number of generations to run each replicate
    %Initial condition: 5 parameters. birth rate, b = (1-m-tau-v), migration rate (m), migration from the regional species pool (v), cladogenesis (tau) and anagenesis (G)
    newspecies = 1;%initialize new species
    R = ones(S,J);%total number of inds in each lake
    m = unifrnd(0.001,0.1,1);
    v = unifrnd(0.0001,0.01,1);%sample from a distribution-pool
    tau = unifrnd(0.0001,0.01,1);%density-indep only tau and K, density-dep tau * (S_T - S_i)/S_i and the strength of density-dep is driven by S_T plus K
    G = unifrnd(0.0001,0.01,1);%gene flow to allow for anagenesis speciation

    %Main loop
    for i = 1:G;%generations
        for j = 1:S*J;%S*J birth-death-migration-speciation loops
            KillHab = unidrnd(S);%pick up lake as a function of its volume
            KillInd = unidrnd(J);%kill ind
            mvb = unifrnd(0,1);%which process: local birth (b), migration (m), pool (m) or cladogenesis (tau)
            if mvb <= m;MigrantHab = unifrnd(0,1);%migration m 
               for k = 1:S;
                   if D(KillHab,k) >= MigrantHab;%check migration matrix M and HM (link D with M), from which lake is coming the migrant?
                      MigrantInd = unidrnd(J);
                      R(KillHab,KillInd) = R(k,MigrantInd);
                   end
                   break
               end
            elseif mvb > mr & mvb <= mr+vr;newspecies = newspecies + 1;%spontaneous == change to tau, initialization of protracted speciation
               R(KillHab,KillInd) = newspecies;                                       
            else
               BirthLocal = unidrnd(J);
               if BirthLocal ~= KillInd;
                  R(KillHab,KillInd) = R(KillHab,BirthLocal);%local birth
               end
               %add probability regional pool
               %Rhone basin new species from the pool: lake number 8 (the closest by distance) or calculate the probability for all the lakes and random chosen lake.
            end
            %check for anagenesis speciation according to G matrix. if some population isolated more generations than G, then new species by anagenesis in lake i
        end%SJ
   end%G
  
%4. Counting species per site
   SpeciesPerSite = zeros(1,4);Beta = zeros(1,2);counti = 0;
   for i = 1:S;
      AR = sort(R(i,:));extantspecies = [ find(AR(1:end-1) ~= AR(2:end)) length(AR) ];richness = AR(extantspecies);
      abu = diff([0 extantspecies]);abu = sort(abu,'descend');
      extantR = length(richness);
   end
   SpeciesPerSite = unique(SpeciesPerSite,'rows');


   %Alpha richness (track SpeciesPerSite)
   count = 0;
   for k = 1:S;
       Alpharichness = find(SpeciesPerSite(:,1) == k);
       count = count + 1;
       Richness(count,1) = length(Alpharichness);
   end

   %Beta and gamma richness (cladogenesis and anagenesis)


%5. Absolute difference (no error control)
   Absdiff = abs(A1(:,1)/(sum(A1)) - Richness/(sum(Richness)));
   zeroall = find(Absdiff == 0);
   Absdiff(Absdiff == 0) = [];
   MLEoctave = log(Absdiff);
   MLET = sum(MLEoctave);

%Absolute difference error control: Data A1; Simu: Richness
count3 = 0;
for u = 1:302;
    if A1(u,1) >= Richness(u,1);
       A3 = 1/((Richness(u,1)+1)*(2 - exp(-1)));
       B3 = exp(-(A1(u,1) - Richness(u,1))/(Richness(u,1)+1));
       count3 = count3 + 1;
       V3(count3,1) = log(A3*B3);
    else
       A3 = 1/((Richness(u,1)+1)*(2 - exp(-1)));
       B3 = exp(-(Richness(u,1) - A1(u,1))/(Richness(u,1)+1));
       count3 = count3 + 1;
       V3(count3,1) = log(A3*B3);
    end
end
TInf = find(V3(:,1) == -Inf);
V3(V3 == -Inf) = [];
VT3 = sum(V3);


%6. plot (testing prototype)
%v=1:302;v=v';
%New = zeros(302,2);
%New(1:302,1) = A1;
%New(1:302,2) = Richness;
%New1 = sort(New,1,'descend')
%plot(v,New1,'o')
%hold on
%plot(v,New1)
%pause

%7. Output file with the pairwise and geographic distance (Distance.m)
fid = fopen('AlphaDataModelParameterMLE.txt','a');fprintf(fid,'%3f %3f %3f %3f %3f %3f %3f\n',v,m,tau,a,1-v-m,G,MLET,length(TInf),VT3);fclose(fid);
end%rr


