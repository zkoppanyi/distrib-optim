% distributed averaging
% Single and multstep gradient methods for averaging
% Papers:
% https://link.springer.com/content/pdf/10.1007/978-0-85729-033-5_4.pdf, Page 128
% http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1307319

 
clear variables; clc;

load('problem')


x0_mat = (ycoors+normrnd(0, 10, 2, size(coors,1)));
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

agents = [];
cont_agents = [size(x0, 1)/2, size(x0, 1)/2-1];

x0i = x0(1:n_params);
d0 = zeros(size(x0i, 1), 1);
for i = 1 : n_agents
        
        negh_i = find(Au(i,:) == 1);
        n_neighb = length(negh_i);
        meas_idx = find(meas(:,1) == i);        
        meas_curr = meas(meas_idx, :);                      

        idxx = ([meas_curr(1,1); meas_curr(:,2)]-1)*2+1;
        idxy = idxx+1;
        x0l_idx = sort([idxx; idxy]);
        x0l = x0i * nan;
        x0l(x0l_idx) = x0i(x0l_idx);       
        n = length(x0l);
        
        r = @(x1, y1, x2, y2) (x1-x2).^2 + (y1-y2).^2;
        r2 = r(x0l(1), x0l(2), x0l(3:2:n), x0l(4:2:n));
        l = meas_curr(:, end);       

        agents{i}.x0l = [x0l, x0l];
        agents{i}.x0l_idx = x0l_idx;
        agents{i}.l = l;
        agents{i}.negh_i = negh_i;
        agents{i}.meas = meas_curr;
end

% init dual variables
for i = 1 : n_agents
    
    negh_i = agents{i}.negh_i;
    n_neighb = length(negh_i);
        
    % Exchange dual variables
    d = zeros(n_params, n_neighb);
    for k2 = 1 : n_neighb 
        didx = intersect(agents{i}.x0l_idx, agents{k2}.x0l_idx);
        d(didx, k2) = d0(didx);  
    end
    agents{i}.d = d;
end

max_iter = 10;
for k1 = 2 : max_iter
    fprintf('iter: %i / %i\n', k1, max_iter )
    for i = 1 : n_agents
        
        negh_i = agents{i}.negh_i;
        n_neighb = length(negh_i);
        meas_idx = find(meas(:,1) == i);        
        meas_curr = meas(meas_idx, :);  
        
        x0l = agents{i}.x0l(:, end);  
        x0l_idx = find(~isnan(x0l));
        x0l2 = x0l(x0l_idx);
        n = length(x0l2);        
        l = agents{i}.l;        
       
        % Exchange dual variables
        d = agents{i}.d;
        dd = zeros(n_params, n_neighb);
        for k2 = 1 : n_neighb
            ni = agents{i}.negh_i(k2);
            d2 = agents{ni}.d;
            idx_curr = find( agents{ni}.negh_i == i );
            if isempty(idx_curr)
                disp('assert error');
                return;
            end
            dd(:,k2) = d(:, k2) - d2(:, idx_curr);  
            
            
        end
        
        sumd = sum(dd,2);
        
        % Solve the dual function
         n = length(x0l);
         x0lo = x0l(1:(n-4));
         idxo = agents{i}.x0l_idx;
         ci = [(cont_agents-1)*2+1, (cont_agents-1)*2+2];
         idxo = setdiff(idxo, ci);        
        
        opts = optimoptions(@fminunc,'Display','none');
        x0le = fminunc(@(x) norm(dual_obj_fn(x, x0i, meas_curr, cont_agents)).^2 + x(idxo)'*sumd(idxo), x0lo, opts);
        x0le = [x0le; x0i((n-3):n)];            
        agents{i}.x0l = [agents{i}.x0l, x0le];               

    end
    
    
    % Ok, let's solve the dual problem and update dual vars
    for i = 1 : n_agents
        negh_i = agents{i}.negh_i;
        n_neighb = length(negh_i);
        
        x0l = agents{i}.x0l(:, end);       
        d = agents{i}.d;
        
        dx0 = zeros(length(d), n_neighb);
        for k2  = 1 : n_neighb 
            x0l2 = agents{negh_i(k2)}.x0l(:, end);
            dx0l2 = x0l - x0l2;
            idx = ~isnan(dx0l2);
            dx0(idx, k2) = dx0l2(idx);
            %dx0l(idx) = dx0l(idx) + 1;
        end
               
        alpha = 0.001;
        de = d + alpha * dx0;
        agents{i}.d = de;        
    end    
    
    
end

x0n = [];
iters_dx = [];
for i = 1 : n_agents
    x0l = agents{i}.x0l(:, end); 
    k = (i-1)*2+1;
    x0n = [x0n; x0l(k), x0l(k+1)];
end

c2 = coors';
iters_dx = agents{1}.x0l - c2(:);

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

