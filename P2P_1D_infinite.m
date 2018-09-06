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
R = 0.1; % effort cost- default is 0.0001
Q = [1 0 0 -1
    0 0 0 0
    0 0 0 0
    -1 0 0 1]*0.1;

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

% [P,eig,L] = dare(Ad,Bd,Q,R);

n = 3000;
P = zeros(order,order,n);
P(:,:,1) = rand(order);
for i = 2:n
    P(:,:,i) = Ad'*P(:,:,i-1)*Ad - (Ad'*P(:,:,i-1)*Bd)*inv(R + Bd'*P(:,:,i-1)*Bd)*(Bd'*P(:,:,i-1)*Ad) + Q;
end
L = inv(R + Bd'*P(:,:,i)*Bd)*(Bd'*P(:,:,i)*Ad);

for i = 2:nstep
    u(:,i) = -L*x(:,i-1);
    x(:,i) = Ad*x(:,i-1) + Bd*u(:,i);
%     actual = C*x(:,i-1);
%     pred = C*xhat(:,i-1);
%     xhat(:,i) = (Ad*xhat(:,i-1) + Bd*u(:,i)) + K*(actual - pred);
%     if i > 2000
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