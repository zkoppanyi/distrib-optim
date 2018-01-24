% distributed averaging
% Single and multstep gradient methods for averaging
% paper: [1] Accelerated Gradient Methods for network Optimizaiton,
% https://arxiv.org/pdf/1211.2132.pdf, Page 14

clear variables; close all; clc;

load('problem')
%xv = rand(n_agents, 1) * 10;
xv = normrnd(15,5,n_agents,1);
mean(xv)

iters = xv;

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
W = L; % weight matrix as graph's laplacian

% 1. Using direct formula
% Au = abs(A);
% W = zeros(n_agents, n_agents);
% alpha = 0.1;
% for i = 1 : n_agents
%     for j = 1 : n_agents
%         if i == j
%             d = sum(Au(i,:));
%             W(i,j) = 1 - d*alpha;
%             continue;
%         end
%         if Au(i,j) == 1
%             W(i, j) = alpha;
%         else
%             W(i,j) = 0;
%         end;
%     end;
% end;

%check conditions
W*ones(size(W,1), 1) 
%== eye(size(W,1), 1)

[S U D] = svd(W-ones(size(W,1), size(W,2))/size(W,1));
max(diag(U)) % < 1

x_prev = xv;
alpha = 0.1;
beta = 0.1;
for k = 1 : 100
    x_next = ((1+beta)*eye(size(W,1), size(W,2)) - alpha*W)*xv - beta*x_prev;
    %x_next = xv - alpha*W*xv;
    
    x_prev = xv;
    xv = x_next;
    
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
xv

plot(iters');
