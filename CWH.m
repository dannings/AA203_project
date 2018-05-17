%% CWH dynamics
% xdot = A*x+B*u, A is already included in Phi
B = [0 0;0 0;1 0;0 1]; % 4-by-2 matrix
%% 2PBVP
%% Setup parameters
mu_E = 3.986004418e5; % gravitational parameter
R_E = 6378; %km
r_p = 6578; r_a = 6878;
a = (r_p+r_a)/2; % semi-major axis
n = sqrt(mu_E/(a)^3);
%% Iterative Search
T = 100:100:2*pi/n; % T = tf-t0, relaxed problem 
x_0 = [100; 300; 0; 0]; % 200 km above the target and zero velocity
x_f = [0.2; 0.7; 0; 0]; % some place near the taget and zero veloctiy
ind_best = 0; % indices of optimal transfer time
dV = zeros(length(T),4); % stores burn delta-Vs
dV_max = 0.5; % km/s
best_cost = Inf;
for i = 1:length(T)
    tau = [0; T(i)];
    Phi_v = cwstm_2d_control(B,n,T(i),tau);
    Phi = cwstm_2d(n, T(i));
    dV(i,:) =Phi_v\(x_f-Phi*x_0);
    dV_1 = dV(i,1:2); dV_2 = dV(i,3:4);
    if norm(dV_1)<=dV_max && norm(dV_2)<=dV_max
        if norm(dV_1)+norm(dV_2)<best_cost
            best_cost = norm(dV_1)+norm(dV_2);
            ind_best = i;
        end
    end
end
tau_optimal = [0; T(ind_best)]; % The optimal time for impulsive manuever
state = zeros(ind_best,4);
state(1,:) = x_0;
% State after initial burn
for i = 2:ind_best-1
    Phi_v_best = cwstm_2d_control(B,n,T(i),tau_optimal(1));
    Phi = cwstm_2d(n,T(i));
    state(i,:) = Phi*x_0+ Phi_v_best*dV(ind_best,1:2)';
end
%State after final burn
Phi_v_best = cwstm_2d_control(B,n,T(ind_best),tau_optimal);
Phi = cwstm_2d(n,T(ind_best));
state(end,:) = Phi*x_0+ Phi_v_best*dV(ind_best,1:4)';
figure(1);
hold on
plot(state(:,2),state(:,1))    
axis tight
grid on
box on
%% Help Function
function Phi =  cwstm_2d(n,dt)
% n: mean motion 
% dt: time counted from soem specific time
%Blocks of matrix Phi
Phi_rr =[4-3*cos(n*dt) 0; 
         6*(sin(n*dt)-n*dt) 1];
Phi_rv = [(1/n)*sin(n*dt) (2/n)*(1-cos(n*dt)); 
          (2/n)*(cos(n*dt)-1) (4/n)*sin(n*dt)-(3/n)*n*dt ];
Phi_vr = [3*n*sin(n*dt) 0;
          6*n*(cos(n*dt)-1) 0];
Phi_vv = [cos(n*dt) 2*sin(n*dt);
          -2*sin(n*dt) 4*cos(n*dt)-3];
Phi = [Phi_rr Phi_rv;
       Phi_vr Phi_vv];
end
function Phi_v = cwstm_2d_control(B,n,t,tau)
% tau is a vector containing times where impulsive velocities are given
Phi_v = zeros(size(B,1),2*length(tau));
for i = 1:length(tau)
    dt = t - tau(i);
    Phi = cwstm_2d(n, dt);
    Phi_v(:,2*i-1:2*i) = Phi*B;
end
end