clear all; clc;

n_agents = 10;
n_neighbors = 5;

if n_agents < n_neighbors 
    disp('incorrect numbers of agents and neighbors.');
    return;
end;

A = zeros(n_agents, n_agents); % incidence matrix
for i = 1 : n_agents
    
    k = 0;
    while (k ~= n_neighbors)
        n = floor(rand * n_agents)+1;
        if n == i,
            continue;
        end;
        %if A(i,n) == 0
            A(i,n) = 1;
            A(n,i) = 1;
            k = k + 1;
        %end;
    end;
    
end;

% is the matrix symmetric?
if (norm(abs(A')-abs(A)) ~=0)
    disp('Connection graph is not symmetric!')
end;

G = graph(abs(A),'OmitSelfLoops')
plot(G);

save('problem')