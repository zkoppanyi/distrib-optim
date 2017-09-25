% Single and multstep gradient methods for quadratic problem 
% paper: [1] Accelerated Gradient Methods for network Optimizaiton,
% https://arxiv.org/pdf/1211.2132.pdf, Page 14

clear all; clc;

load('problem')

f = @(x, xv) (x-xv).^2;
%g = @(x) 1/2*f(x)'*f(x);

%solution
%sol = fminunc(f, 0)


% gradintes
%syms sx1 sx2 sx3
%diff(g([sx1; sx2; sx3]), sx1)

xv = ones(n_agents, 1) + rand(n_agents, 1);
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
W = L; % weight metrix as graph's laplacian

%check conditions
W*ones(size(W,1), 1) 
%== eye(size(W,1), 1)

[S U D] = svd(W-ones(size(W,1), size(W,2))/size(W,1));
max(diag(U)) % < 1

x_prev = xv;
alpha = 0.01;
beta = 0.1;
for k = 1 : 130
    %%x_next = ((1+beta)*eye(size(W,1), size(W,2)) - alpha*W)*xv - beta*x_prev;
    x_next = xv - alpha * W * (2*xv - 6);
    %x_next = xv - alpha * W * xv;
    
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