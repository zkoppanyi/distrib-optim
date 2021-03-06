% distributed averaging
% standard consesus on network-wide mean
% paper: https://pdfs.semanticscholar.org/18ad/bed983ada7e23e31637cb9517436b4cdf392.pdf

clear all; clc;

load('problem')
%xv = rand(n_agents, 1) * 10;
xv = normrnd(15,5,n_agents,1);
mean(xv)

iters = xv;

%Create W
% 1. Using direct formula
Au = abs(A);
W = zeros(n_agents, n_agents);
alpha = 0.1;
for i = 1 : n_agents
    for j = 1 : n_agents
        if i == j
            d = sum(Au(i,:));
            W(i,j) = 1 - d*alpha;
            continue;
        end
        if Au(i,j) == 1
            W(i, j) = alpha;
        else
            W(i,j) = 0;
        end;
    end;
end;

% 2. Using Laplacian 
%create graph Laplacian
Au = abs(A);
C = [];
for i = 1 : n_agents
    for j = (i+1) : n_agents
        if Au(i,j) == 1
            col = zeros(n_agents, 1);
            col(i) = 1;
            col(j) = -1;
            C = [C, col];
        end;
    end;
end;
L = C*C';

[~, S, ~] = svd(L);
%diagS=S(S>1e-5);
diagS = diag(S);
alpha = 2 / (min(diagS) + max(diagS)); %this probably wrong
alpha = 0.1;
W = eye(size(L,1), size(L,2)) - alpha*L;

% Metropolis-Hastings
Au = abs(A);
W = zeros(n_agents, n_agents);
alpha = 0.2;
for i = 1 : n_agents
    di = sum(Au(i,:));    
    for j = 1 : n_agents
       if i ~= j
           if Au(i,j) == 1
               dj = sum(Au(j,:));  
               W(i,j) = min(1/di, 1/dj);
           else
               W(i,j) = 0;
           end
       else
           for k = 1 : n_agents
               if Au(i,k) == 1
                   dk = sum(Au(k,:));  
                   W(i,i) = W(i,i) + max(0, 1/di-1/dk);
               end
           end
       end
       
    end
end

%check conditions
W*ones(size(W,1), 1) 
%== eye(size(W,1), 1)

[S U D] = svd(W-ones(size(W,1), size(W,2))/size(W,1));
max(diag(U)) % < 1

for k = 1 : 60
    xv = W*xv;
    iters = [iters, xv];
end;

%iterations
% Au = abs(A);
% for k = 1 : 60    
%     xv2 = xv;
%     for i = 1 : n_agents
%         sxv = W(i,i)*xv2(i);
%         for j = 1 : n_agents
%             if and(Au(i,j) == 1, i~=j)
%                 sxv = sxv + W(i,j)*xv2(j);
%             end;
%         end;
%         xv(i) = sxv ;
%     end;
%     iters = [iters, xv];
% end;
% xv

plot(iters');
