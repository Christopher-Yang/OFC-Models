% Simulation of a point to point reach with different degrees of rotation.
% This script produces 2D position plots of the three angles of rotation
% the user selects. Parameters to adjust before starting: delt, nstep,
% distance, ang, Q, R.

clear all;
rand('seed',1);

% adjust time step 
delt = 0.001; % time step length in secs
nstep = ceil(3/delt); % number of time steps
distance = -8; % horizontal distance of hand from zero

% values for Q and R taken from Qian infinite horizon model
R = 0.00001; % effort cost- default is 0.0001
Q = [1 0 0 -1
    0 0 0 0
    0 0 0 0
    -1 0 0 1]*0.1;

% % parameters for A and B matrices
% t1 = 0.224;
% t2 = 0.013;
% t3 = 0.004;
% k = 0; %1 to include spring
% b = t1 + t2;
% m = t1*t2;
% r = t3;

% % generate A and B matrices in discrete time formulation
% A = [0 1 0
%     -k/m -b/m 1/m
%     0 0 -1/r];
% A2 = expm(delt*A);
% Ac = blkdiag(A,1);
% Ad = blkdiag(A2,1);
% 
% Bc = [0 0 1/r 0]';
% Bd = delt*Bc;

% Single joint reaching movements:
G = .14;        % Viscous Constant: Ns/m
I = .1;         % Inertia Kgm2
tau = 0.066;    % Muscle time constant, s

A = [0 1 0;0 -G/I 1/I;0 0 -1/tau];
A2 = expm(delt*A);
Ad = blkdiag(A2,1);

B = [0;0;1/tau;0];
Bd = delt*B;

order = size(Ad,1); % order of the system

x = zeros(order,nstep);
x(1,1) = distance; % initialize state variables; x position
xhat = x;
u = zeros(size(Bd,2),nstep); % movement commands

P = zeros(order,order,nstep);
L = zeros(nstep,order);
P(:,:,end) = Q;
for i = 2:nstep
    P(:,:,nstep-i+1) = Ad'*P(:,:,nstep-i+2)*Ad - (Ad'*P(:,:,nstep-i+2)*Bd)*inv(R + Bd'*P(:,:,nstep-i+2)*Bd)*(Bd'*P(:,:,nstep-i+2)*Ad);
    L(nstep-i+1,:) = inv(R + Bd'*P(:,:,nstep-i+2)*Bd)*(Bd'*P(:,:,nstep-i+2)*Ad);
end

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