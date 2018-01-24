% Consensus subgradient method for UWB network optimization
% paper: [1] Subgradient Methods and Consensus Algorithms for Solving Convex Optimization Problems
% http://www.diva-portal.org/smash/get/diva2:496654/FULLTEXT01.pdf, Page 14

clear all; clc;

load('problem')

x0_mat = (ycoors+normrnd(0, 30, 2, size(coors,1)));
%x0_mat = coors';
%fix the last 2 points
x0_mat(1:2, size(x0_mat, 2)) = ycoors(1:2, size(x0_mat, 2));
x0_mat(1:2, size(x0_mat, 2)-1) = ycoors(1:2, size(x0_mat, 2)-1); 
x0r = x0_mat(:);
x0 = x0r;
x0i = x0(1:1:(size(x0_mat, 2)*2-4));
    
iters = [];

% 1. Using direct formula
Au = abs(A);
W1 = zeros(n_agents, n_agents);
alpha = 0.1;
for i = 1 : n_agents
    for j = 1 : n_agents
        if i == j
            d = sum(Au(i,:));
            W1(i,j) = 1 - d*alpha;
            continue;
        end
        if Au(i,j) == 1
            W1(i, j) = alpha;
        else
            W1(i,j) = 0;
        end
    end
end

% 2. Using Laplacian 
% Au = abs(A);
% C = [];
% for i = 1 : n_agents
%     for j = (i+1) : n_agents
%         if Au(i,j) == 1
%             col = zeros(n_agents, 1);
%             col(i) = 1;
%             col(j) = -1;
%             C = [C, col];
%         end
%     end
% end
% L = C*C';
% 
% [~, S, ~] = svd(L);
% %diagS=S(S>1e-5);
% diagS = diag(S);
% %alpha = 2 / (min(diagS) + max(diagS)); %this probably wrong
% alpha = 0.01;
% W = eye(size(L,1), size(L,2)) - alpha*L;

% Metropolis-Hastings
Au = abs(A);
W = zeros(n_agents, n_agents);
alpha = 0.2;
for i = 1 : n_agents
    di = sum(Au(i,:));    
    for j = 1 : n_agents
       if i ~= j
           if Au(i,j) == 1,
               dj = sum(Au(j,:));  
               W(i,j) = min(1/di, 1/dj);
           else
               W(i,j) = 0;
           end;
       else
           for k = 1 : n_agents
               if Au(i,k) == 1,
                   dk = sum(Au(k,:));  
                   W(i,i) = W(i,i) + max(0, 1/di-1/dk);
               end
           end
       end
       
    end
end


%% Check conditions
W*ones(size(W,1), 1);
%== eye(size(W,1), 1)

[S U D] = svd(W-ones(size(W,1), size(W,2))/size(W,1));
max(diag(U)); % < 1


%% Iterations
iters_dx = [];
iters_coors = [];

n_params = size(x0, 1);
%n_params = size(x0, 1);
x0i = [];
x0i= [x0i, x0(1:n_params)];
x0i= [x0i, x0(1:n_params)];
for k1 = 2 : 300

    % Build Jacobian
    u = [];
    for i = 1 : n_agents
        
        idx = find(Au(i,:) == 1);
        n_neighb = length(idx);
        meas_idx = find(meas(:,1) == i);        
        meas_curr = meas(meas_idx, :);               
        
        x0l = x0i(:, k1);
        %for g = 1 : 5,
                        
            J = zeros(size(meas_curr,1), 2*n_agents );
            r = zeros(size(meas_curr,1), 1);
            for k2 = 1 : size(meas_curr,1)

                idx_i = meas_curr(k2,1);
                idx_j = meas_curr(k2,2);
                from_i = (idx_i-1)*2 + 1;
                from_j = (idx_j-1)*2 + 1;

                if Au(idx_i, idx_j) == 1
                    r0 = (x0l(from_i)-x0l(from_j)).^2 + (x0l(from_i+1)-x0l(from_j+1)).^2;            
                    J(k2, from_i)   = (2*x0l(from_j) - 2*x0l(from_i));
                    J(k2, from_i+1) = (2*x0l(from_j+1) - 2*x0l(from_i+1));                       
                    J(k2, from_j)   = (2*x0l(from_i) - 2*x0l(from_j));
                    J(k2, from_j+1) = (2*x0l(from_i+1) - 2*x0l(from_j+1));                    
                    r(k2) = meas_curr(k2, 3)^2-r0;
                end
            end    

            % Removing controls
            from_i = (size(x0_mat, 2)-1)*2 + 1; % last
            J(:, from_i:(from_i+1)) = [];
            from_i = (size(x0_mat, 2)-2)*2 + 1; % last but one
            J(:, from_i:(from_i+1)) = [];

            alpha = 1e-6;
            %dx = alpha * J'*r;
            dx = pinv(J'*J)*J'*r;
            x0l(1:(n_params-4)) = x0l(1:(n_params-4)) - dx;  
        %end
        
        u = [u, x0l];
    end
    
    u = u';   
    u0 = u';
    prev_u = Inf;
    ci = 0;
    %while max(std(u0, [], 2)) > 10e-3
    while max(max(prev_u-u)) > 10e-9
       prev_u = u;
       u = W*u;
       ci = ci+1;
    end
    
    u = u';
    x0i(:, k1+1) = x0i(:, k1);
    x0i(:, k1+1) = u(:,1);
    x0i2 = x0i(:, k1+1);
    x0n = [x0i2(1:2:length(x0i2)), x0i2(2:2:length(x0i2))];
     
    % Errors
    dx = coors-x0n;
    iters_dx = [iters_dx, dx(:)];
    %iters_dx = [iters_dx, u0];
    iters_coors = [iters_coors; x0n'];

end


figure(1); clf; hold on;
plot(iters_dx');
title('Convergence of the parameters');
xlabel('Iterations [-]'); ylabel('Residuals');
set(gca, 'FontSize', 14);

figure(2); clf; hold on;
plot(x0n(:,1), x0n(:,2), 'ro');
plot(coors(:,1), coors(:,2), 'r.', 'MarkerSize', 10);
plot(x0_mat(1,:), x0_mat(2,:), 'b+', 'MarkerSize', 4);
% for k = 1 : 2 : size(iters_coors, 1)
%     lin = iters_coors([k k+1],:)';
%     plot(lin(1,1), lin(1,2), 'b+', 'MarkerSize', 4);
% end
title('Solution');
set(gca, 'FontSize', 14); xlabel('X'); ylabel('Y');
grid on;
axis equal;

    
