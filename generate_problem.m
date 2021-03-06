clear all; clc;

n_agents = 10;
n_neighbors = 6;
type = 2;

if type == 1
    if n_agents < n_neighbors 
        disp('incorrect numbers of agents and neighbors.');
        return;
    end

    A = zeros(n_agents, n_agents); % incidence matrix
    for i = 1 : n_agents

        k = 0;
        while (k ~= n_neighbors)
            n = floor(rand * n_agents)+1;
            if n == i
                continue;
            end
            %if A(i,n) == 0
                A(i,n) = 1;
                A(n,i) = 1;
                k = k + 1;
            %end;
        end
    end

    % is the matrix symmetric?
    if (norm(abs(A')-abs(A)) ~=0)
        disp('Connection graph is not symmetric!')
    end;

    G = graph(abs(A),'OmitSelfLoops')
    plot(G);

end

if type == 2
    step = 100; 
    %sn = floor(sqrt(n_agents));
    sn = n_agents;
    coors = [];
    
    % Regular network
    i = 1;
    for x = 0 : step : sn*step
        for y = 0 : step : sn*step
            coors = [coors; x y];
            i = i+1;
        end
    end
    n_agents = size(coors,1);
    
    % Random network
     coors = rand(n_agents-2, 2)*sn*step;    
     coors = [coors;0 0];
     coors = [coors;sn*step sn*step];
    
    % create connections
    A = zeros(size(coors,1), size(coors,1));
    meas = [];
    for i = 1 : size(coors,1)
        [dist, idx_dist] = sort(sqrt((coors(:,1)-coors(i,1)).^2 + (coors(:,2)-coors(i,2)).^2));
        idx = 2:(n_neighbors+1);
        neigh_idx = idx_dist(idx);
        
        A(i, neigh_idx) = 1;
        A(neigh_idx, i) = 1;
        %meas = [meas; ones(length(idx), 1)*i neigh_idx  dist(idx)];
        meas = [meas; ones(length(idx), 1)*i neigh_idx  dist(idx) + normrnd(0,5)];
    end
    di = logical(eye(size(A)));
    A(di) = 0;

    % is the matrix symmetric?
    if (norm(abs(A')-abs(A)) ~=0)
        disp('Connection graph is not symmetric!')
    end

    % Visualization
    figure(1); clf; hold on;
    %subplot(1, 2, 1); hold on;
    plot(coors(:,1), coors(:,2), 'r.');
    for i = 1 : size(A,1)
        for j = 1 : size(A,1)
            if A(i, j) == 1
                 plot([coors(i,1) coors(j,1)], [coors(i,2) coors(j,2)], 'r-');
            end
            
        end
    end
    xlabel('X'); ylabel('Y');
    axis equal;
    set(gca, 'FontSize', 16);
    
    %subplot(1, 2, 2); hold on;
    %G = graph(abs(A),'OmitSelfLoops');
    %plot(G);
    
    % check the solution
    f = @(x, y, x0, y0, r) (x-x0).^2 + (y-y0).^2 - r.^2;
    norm(f(coors(meas(:,1), 1), coors(meas(:,1), 2), coors(meas(:,2), 1), coors(meas(:,2), 2), meas(:,3))) % has to be zero
    
%     syms x y x0 y0 r
%     diff(f(x, y, x0, y0, r), x)
%     diff(f(x, y, x0, y0, r), y)
%     diff(f(x, y, x0, y0, r), x0)
%     diff(f(x, y, x0, y0, r), y0)
    
    %alpha = 0.2;
    alpha = 3e-7;
    ycoors = coors';
    x0_mat = (ycoors+normrnd(0, 2, 2, size(coors,1)));
    %x0_mat = coors';
    %fix the last 2 points
    x0_mat(1:2, size(x0_mat, 2)) = ycoors(1:2, size(x0_mat, 2));
    x0_mat(1:2, size(x0_mat, 2)-1) = ycoors(1:2, size(x0_mat, 2)-1);   
    
    figure(1);
    %subplot(1, 2, 1); hold on;
    plot(x0_mat(1,:), x0_mat(2,:), 'r.', 'MarkerSize', 15);
    
    x0 = x0_mat(:);
    x0i = x0(1:1:(size(x0_mat, 2)*2-4));
%     for l = 1 : 500
%         
%         % Build Jacobian
%         J = zeros(size(meas,1), size(n_agents*4,1)-4);
%         r = zeros(size(meas,1), 1);
%         for k = 1 : size(meas,1)
% 
%             idx_i = meas(k,1);
%             idx_j = meas(k,2);
%             from_i = (idx_i-1)*2 + 1;
%             from_j = (idx_j-1)*2 + 1;
% 
% %             r0 =  sqrt((x0(from_i)-x0(from_j)).^2 + (x0(from_i+1)-x0(from_j+1)).^2);            
% %             J(k, from_i) = (x0(from_j) - x0(from_i)) / r0;
% %             J(k, from_i+1) = (x0(from_j+1) - x0(from_i+1)) / r0;                       
% %             J(k, from_j) = (x0(from_i) - x0(from_j)) / r0;
% %             J(k, from_j+1) = (x0(from_i+1) - x0(from_j+1)) / r0;                    
% %             r(k) = meas(k, 3)-r0;
% 
%             r0 =  (x0(from_i)-x0(from_j)).^2 + (x0(from_i+1)-x0(from_j+1)).^2;            
%             J(k, from_i)   = (2*x0(from_j) - 2*x0(from_i));
%             J(k, from_i+1) = (2*x0(from_j+1) - 2*x0(from_i+1));                       
%             J(k, from_j)   = (2*x0(from_i) - 2*x0(from_j));
%             J(k, from_j+1) = (2*x0(from_i+1) - 2*x0(from_j+1));                    
%             r(k) = meas(k, 3)^2-r0;
% 
%         end
%         from_i = (size(x0_mat, 2)-1)*2 + 1;
%         J(:, from_i:(from_i+1)) = [];
%         from_i = (size(x0_mat, 2)-2)*2 + 1;
%         J(:, from_i:(from_i+1)) = [];
%                 
% %         if (l > 450)
% %            x0i = x0i - inv(J'*J)*J'*r; % gauss-newton
% %         else
% %             x0i = x0i - alpha*J'*r; % gradient
% %         end
%         
%         x0(1:(size(x0_mat, 2)*2-4)) = x0i(1:(size(x0_mat, 2)*2-4));
%         %norm(J'*r)
%         %norm(x0 - ycoors(:))
%         
%         x0n = [x0(1:2:length(x0)), x0(2:2:length(x0))];
%         norm(f(x0n(meas(:,1), 1), x0n(meas(:,1), 2), x0n(meas(:,2), 1), x0n(meas(:,2), 2), meas(:,3)))
%     end
    
    cont_agents = [size(x0, 1)/2, size(x0, 1)/2-1];
    ci = [(cont_agents-1)*2+1, (cont_agents-1)*2+2];        
    %opts = optimoptions(@fminunc,'Display','none');
    %n = length(x0i);
    %x0i0 = x0i(1:(n-4));
    chk_prob.cont_agents = cont_agents;
    chk_prob.x0 = x0;
    x0i = lsqnonlin(@(x) obj_fn(x, chk_prob.x0, meas, chk_prob.cont_agents), x0i);
    chk_prob.res =  norm(obj_fn(x0i, chk_prob.x0, meas, chk_prob.cont_agents));
        
    % rebuild
    n = length(x0);
    x0n1 = [x0i(1:2:length(x0i)), x0i(2:2:length(x0i))];    
    x0i = [x0i; x0((n-3):n)];
    x0n = [x0i(1:2:length(x0i)), x0i(2:2:length(x0i))];    
    figure(1);
    %subplot(1, 2, 1); hold on;
    plot(x0n1(:,1), x0n1(:,2), 'b.');
    plot(x0n1(:,1), x0n1(:,2), 'bo');
    xlima = [min(coors(:,1))-100, max(coors(:,1))+100];
    ylima = [min(coors(:,2))-100, max(coors(:,2))+100];
    xlim(xlima); ylim(ylima);
    grid on;
end
norm(x0n(:) - coors(:))

sol = x0n;
clear x0n;
x0_mat = (sol'+normrnd(0, 25, 2, size(coors,1)));
save('problem')


% is graph biparite

% events = {'edgetonew', 'edgetodiscovered', 'edgetofinished'};
% T = bfsearch(G, 1, events, 'Restart', true);
% partitions = false(1, numnodes(G));
% is_bipart = true;
% 
% for ii=1:size(T, 1)   
%     if T.Event(ii) == 'edgetonew'
%         partitions(T.Edge(ii, 2)) = ~partitions(T.Edge(ii, 1));
%     else
%         if partitions(T.Edge(ii, 1)) == partitions(T.Edge(ii, 2))
%             is_bipart = false;
%             break;
%         end
%     end
% end
% is_bipart