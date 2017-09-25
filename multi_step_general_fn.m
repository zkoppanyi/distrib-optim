% Single and multstep gradient methods for quadratic problem 
% paper: [1] Accelerated Gradient Methods for Networked Optimizaiton,
% https://arxiv.org/pdf/1211.2132.pdf, Page 14

clear all; clc;

load('problem')

x0_mat = (ycoors+normrnd(0, 0.1, 2, size(coors,1)));
%x0_mat = coors';
%fix the last 2 points
x0_mat(1:2, size(x0_mat, 2)) = ycoors(1:2, size(x0_mat, 2));
x0_mat(1:2, size(x0_mat, 2)-1) = ycoors(1:2, size(x0_mat, 2)-1); 
x0r = x0_mat(:);
x0 = x0r;
x0i = x0(1:1:(size(x0_mat, 2)*2-4));
    
iters = [];

% 2. Using Laplacian 
%create graph Laplacian
Au = abs(A);
C = zeros(n_agents*2, n_agents*2);
for i = 1 : n_agents
    for j = (i+1) : n_agents
        if Au(i,j) == 1
            col = zeros(n_agents*2, 1);
            col((i-1)*2+1) = 1;
            col((i-1)*2+2) = 1;
            col((j-1)*2+1) = 1;
            col((j-1)*2+2) = 1;
            C = [C, col];
        end
    end
end
% C(:,end) = []; C(:,end) = []; C(:,end) = []; C(:,end) = [];
% C(end,:) = []; C(end,:) = []; C(end,:) = []; C(end,:) = [];

L = C*C';
W = L; % weight metrix as graph's laplacian

%check conditions
W*ones(size(W,1), 1) 
%== eye(size(W,1), 1)

[S U D] = svd(W-ones(size(W,1), size(W,2))/size(W,1));
max(diag(U)) % < 1

x0i = x0;
x_prev = x0i;
alpha = 0.01;
beta = 0.001;
for k = 1 : 130
    
    % Build Jacobian
    J = zeros(size(meas,1), size(n_agents*4,1));
    r = zeros(size(meas,1), 1);
    for k = 1 : size(meas,1)

        idx_i = meas(k,1);
        idx_j = meas(k,2);
        from_i = (idx_i-1)*2 + 1;
        from_j = (idx_j-1)*2 + 1;

        r0 =  sqrt((x0(from_i)-x0(from_j)).^2 + (x0(from_i+1)-x0(from_j+1)).^2);

        J(k, from_i) = (x0(from_j) - x0(from_i)) / r0;
        J(k, from_i+1) = (x0(from_j+1) - x0(from_i+1)) / r0;           

        J(k, from_j) = (x0(from_i) - x0(from_j)) / r0;
        J(k, from_j+1) = (x0(from_i+1) - x0(from_j+1)) / r0;        

        r(k) = meas(k, 3)-r0;
    end
%     from_i = (size(x0_mat, 2)-1)*2 + 1;
%     J(:, from_i:(from_i+1)) = 0;
%     from_i = (size(x0_mat, 2)-2)*2 + 1;
%     J(:, from_i:(from_i+1)) = 0;
        
    %x_next = ((1+beta)*eye(size(W,1), size(W,2)) - alpha*W)*J'*r - beta*(x_prev - x0i);
    x_next = x0i - alpha * W * J'*r;
    
    x_prev = x0i;
    x0i = x_next;
    x0(1:(size(x0_mat, 2)*2-4)) = x0i(1:(size(x0_mat, 2)*2-4));
     
    iters = [iters, x0(:)-x0r];

    x0n = [x0(1:2:length(x0)), x0(2:2:length(x0))];
    norm(f(x0n(meas(:,1), 1), x0n(meas(:,1), 2), x0n(meas(:,2), 1), x0n(meas(:,2), 2), meas(:,3)))
end

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
x0

figure(1); clf; hold on;
plot(iters');

figure(2); clf; hold on;
plot(x0n(:,1), x0n(:,2), 'ro');
plot(x0_mat(1,:), x0_mat(2,:), 'b*');
plot(coors(:,:), coors(:,1), 'r*');

    
