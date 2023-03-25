clc;
close all;
clear;

% Define system parameters
C1 = 100; % Heat capacity of room 1 (J/K)
C2 = 150; % Heat capacity of room 2 (J/K)
h = 0.05; % Heat transfer coefficient btw room and outside temp (W/K.m^2)
A_s = 20; % Surface area of room 1 (m^2)
A_s2 = 40; % Surface area of room 2 (m^2)
k = 0.7; % Heat transfer coefficient btw the rooms(W/K.m^2)
Tout = 20; % Outside temperature (degC)
T_desired = [100; 150]; % Desired temperature of each room (degC)


%System State Space Representation
A = [-(h*A_s+k)/C1, k/C1;
    k/C2, -(h*A_s2+k)/C2];
B = [1/C1 0; 0, 1/C2];
C = [1 0; 0 1];
D = 0*eye(2);


% Calculate the observability matrix
CO = ctrb(A,B);

% Check if the system is controllable
if rank(CO) == size(A,1)
    disp('System is controllable!');
else
    disp('System is not controllable!');
end

% Calculate the observability matrix
O = obsv(A, C);

% Check the rank of the observability matrix
if rank(O) == size(A,1)
    disp('System is observable');
else
    disp('System is not observable');
end


%Define cost function
Q1 = diag([1,1]);
Q2 = 2*(C'*C);
R1 = eye(2);
R2 = eye(2);

% LQR controller
[K1,S1,P1] = lqr(A,B,Q1,R1);
[K2,S2,P2] = lqr(A,B,Q2,R2);


% Kalman filter
G = eye(2); % Process noise covariance matrix
Q = 1*eye(2); % Measurement noise covariance matrix
R = 100*eye(2); % Error covariance matrix
[L,P] = lqe(A,G,C,Q,R); % Kalman gain

% Define simulation parameters
tspan = [0 3500]; % Simulation time span (s)
initial_state = [20; 20]; % Initial temperature of each room (degC)
u1 = 100; % Heat input to room 1 (W)
u2 = 150; % Heat input to room 2 (W)

% Define the system dynamics as a function of time and state
f = @(t, S) [(-h*A_s/C1)*(S(1)-Tout) + k/C1*S(2) + u1/C1;
k/C2*S(1) + (-h*A_s2/C2)*(S(2)-Tout) + u2/C2];

% Simulate the system using ODE45 solver
[t, S] = ode45(f, tspan, initial_state);


% Plot the results
figure;
plot(t, S(:,1), 'b-', t, S(:,2), 'r-');
xlabel('Time (s)');
ylabel('Temperature (degC)');
legend('Room 1', 'Room 2', 'Location', 'southeast');
title('Simulation of two-room model');



%% Specify disturbance and noise magnitude
% Vd = eye(2); % disturbance covariance
% Vn = eye(2); % noise covariance
% 
% % Build Kalman filter
% [Kf,P,E] = lqe(A,eye(2),C,Vd,Vn); % design Kalman filter
% % alternatively, possible to design using "LQR" code
% Kq = (lqr(A',C',Vd,Vn))';
% 
% %% Augment system with additional inputs
% B_aug = [B eye(2) 0*B]; % [u I*wd 0*wn]
% D_aug = [0 0 0 0 0 0; 0 0 0 0 0 1]; % D matrix passes noise through
% sysC = ss(A,B_aug,C,D_aug); % single-measurement system
% 
% % "true" system w/ full-state output, disturbance, no noise
% sysTruth = ss(A,B_aug,eye(2),zeros(2,size(B_aug,2)));
% sysKF = ss(A-Kf*C,[B Kf],eye(2),0*[B Kf]); % Kalman filter


% dt = .01;
% t_kf = dt:dt:50;
% uDIST = sqrt(Vd)*randn(2,size(t_kf,2)); % random disturbance
% uNOISE = sqrt(Vn)*randn(size(t_kf,2)); % random noise
% u = 0*t;
% u(1/dt) = 20/dt; % positive impulse
% u(15/dt) = -20/dt; % negative impulse
% u_aug = [u; uDIST; uNOISE]; % input w/ disturbance and noise
% [y,t_kf] = lsim(sysC,u_aug,t_kf); % noisy measurement
% [xtrue,t_kf] = lsim(sysTruth,u_aug,t_kf); % true state
% [xhat,t_kf] = lsim(sysKF,[u; y'],t_kf); % state estimate
