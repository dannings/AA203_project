function [c, N, A] = FMTstar(start, goal, obstacles, Nmax, eta, FOR_THE_KIDDOS)
% Please refer to the pseudo code in the paper
clf; hold on;
axis square; rectangle;
plot(start(1), start(2), 'sg')
fill(goal(1,:), goal(2,:), 'g') % plot obstacles
mu = 1;
for i = 1:length(obstacles)
    fill(obstacles{i}(1,:), obstacles{i}(2,:), 'r')
    mu = mu - polyarea(obstacles{i}(1,:), obstacles{i}(2,:)); % remaining area that is free of obstacles
end

%% Sample Free
V = [start, rand(2,Nmax)]; % generate sample points inside the region, include start point, Line 1
for i = 1:length(obstacles)
    V = V(:,arrayfun(@(x,y) ~inpolygon(x, y, obstacles{i}(1,:), obstacles{i}(2,:)), V(1,:), V(2,:))); 
    % check if the points in the generated batch is in V_obstacle, throw
    % them away if they are in obstacles, obtaining V_free
end
N = size(V,2); % number of  free sample points 
%======================Need to be modified¡ý=================================
D = squareform(pdist(V'));  
% Distance between pairs of free sample points, 0 at the diagonal;
% In our case need to be the minimum fuel consumption
%======================Need to be modified¡ü=================================
A = zeros(1, N); % Storing minimum value (y) for each section
plot(V(1,:),V(2,:),'ko');
%%
%======================Need to be modified¡ý=================================
r = eta*sqrt(2*mu/pi*log(N)/N); % This can be changed 
% NN stores the set of samples{u ¡Ê V_free : |u-v|<r} for each point(NN is a
% cell), line 4
NN = arrayfun(@(v) [find(D(v,1:v-1) < r), v + find(D(v,v+1:end) < r)], 1:N, 'uniformoutput', false);
% The D is different in our case (not distance)
%======================Need to be modified¡ü=================================
W = true(1,N);
W(1) = false; % V_unvisited boolean, excluding x_init
H = false(1,N);
H(1) = true; % V_open boolean, only includes x_init
C = inf*ones(1,N);
C(1) = 0;
z = 1; % index, Line 3 in algorithm of paper

while ~inpolygon(V(1,z), V(2,z), goal(1,:), goal(2,:)) % while z not in region V_goal
    H_new = []; % V_open_new
    for x = NN{z}(W(NN{z})) % Line 8 & 9 For each x in X_near (X_near is consist of unvisited 'near' points)
        Y_near = NN{x}(H(NN{x})); % Line 10, 11 and 12 ,'near'points in open set around this x 
        %======================Need to be modified¡ý=================================
        [c_min, y_idx] = min(C(Y_near) + D(x, Y_near)); % Here we insert the dynamics instea dof D
        %======================Need to be modified¡ü=================================
        y_min = Y_near(y_idx); % find the "near" point in open set that gives minimal cost, Line 13  
        if check_segment(V(:,x), V(:,y_min), obstacles) % Collision Check, If collision free, Line 14
            A(x) = y_min;
            C(x) = c_min;
            %======================Need to be modified¡ý=================================
            plot(V(1,[x, y_min]), V(2,[x, y_min]), 'ok-')
            % Need to calculate path through CWH equations
            %======================Need to be modified¡ü=================================
            if FOR_THE_KIDDOS
                drawnow;
            end
            H_new(end+1) = x; % add current x to V_open_new
            W(x) = false; % remove x from V_unvisited
        end
    end
    H(H_new) = true; % V_open = V_open combined with V_open_new
    H(z) = false; % Exclude z from V_open, Line 21-22
    z = find(C == min(C(H)));   % z = argmin(y included in V_open){c(y)}, Line 26
end

v = z;
X = V(1,v);
Y = V(2,v);
while v ~= 1
    v = A(v);
    X(end+1) = V(1,v);
    Y(end+1) = V(2,v);
end
plot(X,Y,'ob-','LineWidth',2,'MarkerFaceColor','b');
hold off;
c = C(z);
% cumsum(fliplr(sqrt(diff(X).^2 + diff(Y).^2)))
end
%% Help function
function valid = check_segment(a, b, obstacles) % Check Collision
valid = all(cellfun(@(obs) isempty(polyxpoly([a(1) b(1)], [a(2) b(2)], obs(1,:), obs(2,:))), obstacles));
end