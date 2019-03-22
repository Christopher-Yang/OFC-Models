% Simulation of a 1D point-to-point reach in a finite horizon formulation.
% This script produces position, velocity, and motor command plots for the
% simulation. Parameters to adjust before starting: delt, simTime, Q, R.

clear all;
rand('seed',1);

delt = 0.001; % time step length in secs
simTime = 3; % number of seconds to simulate movements
nstep = ceil(simTime/delt); % number of time steps

% accuracy and effort matrices
R = 0.00001; % effort cost- default is 0.0001
Q = [1 0 0 -1
    0 0 0 0
    0 0 0 0
    -1 0 0 1]*0.1;

% Single joint reaching movements:
G = .14;        % Viscous Constant: Ns/m
I = .1;         % Inertia Kgm2
tau = 0.066;    % Muscle time constant, s

A = [0 1 0;0 -G/I 1/I;0 0 -1/tau]; % build system dynamics matrix
A2 = expm(delt*A); % discretize A
Ad = blkdiag(A2,1); % augment A with target location dynamics

B = [0;0;1/tau;0]; % build input matrix
Bd = delt*B; % discretize B

order = size(Ad,1); % order of the system

x = zeros(order,nstep);
x(1,1) = 0; % hand starting position
x(4,1) = 3;
u = zeros(size(Bd,2),nstep); % movement commands

% calculate control law, L
P = zeros(order,order,nstep);
L = zeros(nstep,order);
P(:,:,end) = Q;
for i = 2:nstep
    P(:,:,nstep-i+1) = Ad'*P(:,:,nstep-i+2)*Ad - (Ad'*P(:,:,nstep-i+2)*Bd)*inv(R + Bd'*P(:,:,nstep-i+2)*Bd)*(Bd'*P(:,:,nstep-i+2)*Ad);
    L(nstep-i+1,:) = inv(R + Bd'*P(:,:,nstep-i+2)*Bd)*(Bd'*P(:,:,nstep-i+2)*Ad);
end

% simulate arm movements
for i = 2:nstep
    u(:,i) = -L(i-1,:)*x(:,i-1);
    x(:,i) = Ad*x(:,i-1) + Bd*u(:,i);
%     if i > 50000
%         x(4,i) = -2;
%     end
end

figure
subplot(1,3,1)
plot(0:delt:delt*nstep-delt,x(1,:))
title('Position')
xlabel('Time (s)')
subplot(1,3,2)
plot(0:delt:delt*nstep-delt,x(2,:))
title('Velocity')
xlabel('Time (s)')
subplot(1,3,3)
plot(0:delt:delt*nstep-delt,u)
title('Motor Commands')
xlabel('Time (s)')