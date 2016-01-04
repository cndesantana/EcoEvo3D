%CJ Melian@EAWAG, 3d-dendritic network
DendriticNetworkED;Dendritic;
HM = zeros(8,8);%hydrological Matrix (3d)
M = zeros(8,8);%Migration Matrix
n = 0.001;%Upstream cost
Q = B + B';

%from i to j
for i = 1:7;
     for j = i+1:8;k = zeros(1,1);
        [dist p] = dijkstra(Q,i,j);
        d = zeros(1,length(p) - 1);
        %fill hydrological distance matrix
        for k = 1:length(p) - 1;
            if A(p(1,k),2) > A(p(1,k+1),2);
               A(p(1,k),2);
               A(p(1,k+1),2);
               d(1,k) = Q(p(1,k),p(1,k+1));
            else
                A(p(1,k),2);
                A(p(1,k+1),2);
               d(1,k) = Q(p(1,k),p(1,k+1))^((1 + n*(A(p(1,k+1),2) - A(p(1,k),2)))/1);
            end
        end
        TD = sum(d);
        HM(i,j) = TD;
     end
end
%from j to i
for i = 1:7;
     for j = i+1:8;k = zeros(1,1);
        [dist p] = dijkstra(Q,i,j);
        d = zeros(1,length(p) - 1);
        %fill hydrological distance matrix
        for k = 1:length(p) - 1;
            if A(p(1,k),2) < A(p(1,k+1),2);
               A(p(1,k),2);
               A(p(1,k+1),2);
               d(1,k) = Q(p(1,k),p(1,k+1));
            else
                A(p(1,k),2);
                A(p(1,k+1),2);
               d(1,k) = Q(p(1,k),p(1,k+1))^((1 + n*(A(p(1,k),2) - A(p(1,k+1),2)))/1);
            end
        end
        TD = sum(d);
        HM(j,i) = TD;
        
     end
end
%import volume each lake
%Migration matrix
count = 0;
for j = 1:8;
    M(j,:) = A(1:8,3)
    count = count + 1;
    M(count,count) = 0;
end

