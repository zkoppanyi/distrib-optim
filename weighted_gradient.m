% Weighted Gradient Method for UWB network optimization
% paper: [1] Accelerated Gradient Methods for Networked Optimizaiton,
% https://arxiv.org/pdf/1211.2132.pdf, Page 14

clear all; clc;

load('problem')

x0_mat = (sol'+normrnd(0, 40, 2, size(coors,1)));

%x0_mat = coors';
%fix the last 2 points
x0_mat(1:2, size(x0_mat, 2)) = ycoors(1:2, size(x0_mat, 2));
x0_mat(1:2, size(x0_mat, 2)-1) = ycoors(1:2, size(x0_mat, 2)-1); 
x0r = x0_mat(:);
x0 = x0r;
x0i = x0(1:1:(size(x0_mat, 2)*2-4));

%% Using Laplacian
% Create graph Laplacian
Au = abs(A);
C = zeros(n_agents*2, n_agents*2);
for i = 1 : n_agents
    for j = (i+1) : n_agents
        if Au(i,j) == 1
            col = zeros(n_agents*2, 1);
            col((i-1)*2+1) = 1;
            col((i-1)*2+2) = 1;
            col((j-1)*2+1) = -1;
            col((j-1)*2+2) = -1;
            C = [C, col];
        end
    end
end

L = C*C';

n = 2*max(diag(L)); % max deg.
epsilon = 1/n/2;
epsilon = epsilon / 2;

%n = sum(diag(L)); 
%epsilon = 1/n;

%epsilon = 1e-9;
W = eye(size(L,1)) - epsilon*L; % weight matrix as graph's laplacian

Wt = W;
for k = 1 : 100
    Wt = Wt * Wt;
end


%% 4. Minimization
% % from https://pdfs.semanticscholar.org/18ad/bed983ada7e23e31637cb9517436b4cdf392.pdf, Page 69.
% % Need the Laplacian
% n = size(C,1);
% m = size(C,2);
% %obj_fn = @(w) norm(eye(n) - C*diag(w)*C' - ones(n)/n);
% obj_fn = @(w) norm(eye(n) - C*diag(w)*C');
% w0 = ones(m,1);
% w = fminunc(obj_fn, w0);
% W = eye(n) - C*diag(w)*C';

%% Check conncetion matrix

%check conditions
W*ones(size(W,1), 1);
%== eye(size(W,1), 1)

%v = eig(W-ones(size(W,1), size(W,2))/size(W,1));
v = eig(W-1/size(W,1));
max(v) % < 1

% is symmetric 
sum(sum(W' - W)) == 0;

%% Iterations
iters_dx = [];
iters_coors = [];

n_params = size(x0, 1)-4;
%n_params = size(x0, 1);
x0i = [];
x0i = [x0i, x0(1:n_params)];
x0i = [x0i, x0(1:n_params)];
dit = [];
for k1 = 2 : 600

    % Build Jacobian
    J = zeros(size(meas,1), 2*n_agents );
    r = zeros(size(meas,1), 1);
    for k2 = 1 : size(meas,1)

        idx_i = meas(k2,1);
        idx_j = meas(k2,2);
        from_i = (idx_i-1)*2 + 1;
        from_j = (idx_j-1)*2 + 1;

%         r0 =  sqrt((x0(from_i)-x0(from_j)).^2 + (x0(from_i+1)-x0(from_j+1)).^2);
%         J(k2, from_i) = (x0(from_j) - x0(from_i)) / r0;
%         J(k2, from_i+1) = (x0(from_j+1) - x0(from_i+1)) / r0;         
%         J(k2, from_j) = (x0(from_i) - x0(from_j)) / r0;
%         J(k2, from_j+1) = (x0(from_i+1) - x0(from_j+1)) / r0;        
%         r(k2) = meas(k2, 3)-r0;
        
        r0 = (x0(from_i)-x0(from_j)).^2 + (x0(from_i+1)-x0(from_j+1)).^2;            
        J(k2, from_i)   = (2*x0(from_j) - 2*x0(from_i));
        J(k2, from_i+1) = (2*x0(from_j+1) - 2*x0(from_i+1));                       
        J(k2, from_j)   = (2*x0(from_i) - 2*x0(from_j));
        J(k2, from_j+1) = (2*x0(from_i+1) - 2*x0(from_j+1));                    
        r(k2) = meas(k2, 3)^2-r0;

    end    
    
    % Removing controls
    W2 = W;
    from_i = (size(x0_mat, 2)-1)*2 + 1; % last
    J(:, from_i:(from_i+1)) = [];
    W2(:, from_i:(from_i+1)) = [];
    W2(from_i:(from_i+1), :) = [];

    from_i = (size(x0_mat, 2)-2)*2 + 1; % last but one
    J(:, from_i:(from_i+1)) = [];
    W2(:, from_i:(from_i+1)) = [];
    W2(from_i:(from_i+1), :) = [];
        
    % check solvability 
    if min(eig(J'*J)) <= 0
        fprintf('Iteration #%i\n', k1);
        fprintf('Not solvable! Hessian not positive definite\n');
        %fprintf('Min deg of verteces: %i\n', min(diag(L)));
        break;
    end
    
    
    % Accelerated version
    if k1 == 2
        %alpha = 0.003; beta = 0.1;
        H = J'*J;
        lams = eig(W2*H);
        lams = lams(lams~=0); % remove zero eigens
        min_lam = min(lams);
        max_lam = max(lams);
        alpha = (2 / (sqrt(max_lam) + sqrt(min_lam)) )^2 /10;
        beta = ( (sqrt(max_lam) - sqrt(min_lam)) / (sqrt(max_lam) + sqrt(min_lam)) )^2 / 10; %optimal parameters
    end
    
    alpha = 8e-7;
    beta = 1e-1;
    x0i(:, k1+1) = x0i(:, k1) - alpha * W2 * J'*r + beta*(x0i(:, k1) - x0i(:, k1-1));

    % Simple version
    %alpha = 0.003;
    %x0i(k1+1) = x0i(k1) - alpha * W * J'*r;
    
    x0(1:n_params) = x0i(:, k1+1);    
    x0n = [x0(1:2:length(x0)), x0(2:2:length(x0))];
    
    % Errors
    dx = x0n - coors;
    iters_dx = [iters_dx, dx(:)];
    iters_coors = [iters_coors; x0n'];
    norm(f(x0n(meas(:,1), 1), x0n(meas(:,1), 2), x0n(meas(:,2), 1), x0n(meas(:,2), 2), meas(:,3)));
        
    r = obj_fn(x0, chk_prob.x0, meas, chk_prob.cont_agents);    
    dit =[dit, ((norm(r) / chk_prob.res) - 1)*100]; 
    %dit =[dit, n2]; 
    %if ( abs(n2-n1)/n1*100 < 0.5 )
    %    %break;
    %end
end


figure(1); clf; hold on;
yyaxis right
plot(dit', 'LineWidth', 3);
ylabel('(norm(r)-norm(r^*))/norm(r^*) [%]');
idx = find(abs(diff(dit)) < 0.01);
if ~isempty(idx)
    plot([idx(1) idx(1)], [0 max(dit)], 'r--', 'LineWidth', 2);
end
dit(idx(1));
txt = sprintf('#%i, %.1f%%', idx(1), dit(idx(1)))
text(idx(1)+10, dit(idx(1))+50, txt,'FontSize',14,'Color','red');
xlabel('Iterations [-]'); 
set(gca, 'FontSize', 14);

yyaxis left
plot(iters_dx', '-');
ylabel('Residuals');
xlabel('Iterations [-]'); 
set(gca, 'FontSize', 14);

figure(2); clf; hold on;
plot(coors(:,1), coors(:,2), 'r.', 'MarkerSize', 15);
for i = 1 : size(A,1)
        for j = 1 : size(A,1)
            if A(i, j) == 1
                 plot([coors(i,1) coors(j,1)], [coors(i,2) coors(j,2)], 'r-');
            end
            
        end
end
plot(x0n(:,1), x0n(:,2), 'b.');
plot(x0n(:,1), x0n(:,2), 'bo');
    
plot(x0_mat(1,:), x0_mat(2,:), 'b+', 'MarkerSize', 4);
% for k = 1 : 2 : size(iters_coors, 1)
%     lin = iters_coors([k k+1],:)';
%     plot(lin(1,1), lin(1,2), 'b+', 'MarkerSize', 4);
% end
%title('Solution');
set(gca, 'FontSize', 14); xlabel('X'); ylabel('Y');
grid on;
axis equal;
xlim(xlima); ylim(ylima);

    
