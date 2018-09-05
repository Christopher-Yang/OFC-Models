% Simulation of a point to point reach with different degrees of rotation.
% This script produces 2D position plots of the three angles of rotation
% the user selects. Parameters to adjust before starting: delt, nstep,
% distance, ang, Q, R.

clear all;
rand('seed',1);

% adjust time step 
delt = 0.00001; % time step length in secs
nstep = ceil(1/delt); % number of time steps
distance = -8; % horizontal distance of hand from zero

% values for Q and R taken from Qian infinite horizon model
R = 0.1; % effort cost- default is 0.0001
Q = [1 0 0 -1
    0 0 0 0
    0 0 0 0
    -1 0 0 1]*0.1;

% parameters for A and B matrices
t1 = 0.224;
t2 = 0.013;
t3 = 0.004;
k = 0; %1 to include spring
b = t1 + t2;
m = t1*t2;
r = t3;

% generate A and B matrices in discrete time formulation
A = [0 1 0
    -k/m -b/m 1/m
    0 0 -1/r];
A2 = expm(delt*A);
Ac = blkdiag(A,1);
Ad = blkdiag(A2,1);

Bc = [0 0 1/r 0]';
Bd = delt*Bc;
C = [1 0 0 0
    0 0 0 1];

order = size(Ad,1); % order of the system

x = zeros(order,nstep);
x(1,1) = distance; % initialize state variables; x position
xhat = x;
u = zeros(size(Bd,2),nstep); % movement commands

% [P,eig,L] = dare(Ad,Bd,Q,R);

n = 180000;
P = zeros(order,order,n);
P(:,:,1) = rand(order);
for i = 2:n
    P(:,:,i) = Ad'*P(:,:,i-1)*Ad - (Ad'*P(:,:,i-1)*Bd)*inv(R + Bd'*P(:,:,i-1)*Bd)*(Bd'*P(:,:,i-1)*Ad) + Q;
end
L = inv(R + Bd'*P(:,:,i)*Bd)*(Bd'*P(:,:,i)*Ad);

sys = ss(Ad,Bd,C,0);
Qn = 0.0001; % process noise
Rn = diag(repmat(0.001,[2 1])); % measurement noise
% [kest, K, P2] = kalman(sys,Qn,Rn,0);

for i = 2:nstep
    u(:,i) = -L*x(:,i-1);
    x(:,i) = Ad*x(:,i-1) + Bd*u(:,i);
%     actual = C*x(:,i-1);
%     pred = C*xhat(:,i-1);
%     xhat(:,i) = (Ad*xhat(:,i-1) + Bd*u(:,i)) + K*(actual - pred);
    if i > 70000
        x(4,i) = -2;
    end
end

plot(x(1,:));