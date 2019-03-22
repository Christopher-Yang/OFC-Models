% Simulation of 1D point to point reach in an inifinite horizon 
% formulation. This script produces position and velocity plots of the 
% simulated reach. Parameters to adjust before starting: delt, simTime, Q, 
% R.

clear all;
rng(1);

delt = 0.001; % time step length in secs
simTime = 1; % number of seconds to simulate movements
nstep = ceil(simTime/delt); % number of time steps

% accuracy and effort matrices
R = 0.0001; % effort cost- default is 0.0001
Q = [1.2 0 0 -1
    0 0.02 0 0
    0 0 0 0
    -1 0 0 1];

% Single joint reaching movements:
G = .14;        % Viscous Constant: Ns/m
I = .1;         % Inertia Kgm2
tau = 0.066;    % Muscle time constant, s

A = [0 1 0;0 -G/I 1/I;0 0 -1/tau]; % system dynamics matrix
A2 = expm(delt*A); % discretize A
Ad = blkdiag(A2,1); % augment matrix to include target 

B = [0;0;1/tau;0]; % input matrix
Bd = delt*B; % discretize B

order = size(Ad,1); % determine the order of the system

% initialize the state vector, x
x = zeros(order,nstep); % starting velocity and acceleration set to 0
x(1,1) = 0; % set starting hand position
x(4,1) = 5; % set ending hand position
u = zeros(size(Bd,2),nstep); % movement commands

% calculate control law, L
n = 3000; % number of times to iterate over P 
P = zeros(order,order,n);
P(:,:,1) = rand(order);
for i = 2:n
    P(:,:,i) = Ad'*P(:,:,i-1)*Ad - (Ad'*P(:,:,i-1)*Bd)*inv(R + Bd'*P(:,:,i-1)*Bd)*(Bd'*P(:,:,i-1)*Ad) + Q;
end
L = inv(R + Bd'*P(:,:,i)*Bd)*(Bd'*P(:,:,i)*Ad); % control law

% [P,eig,L] = dare(Ad,Bd,Q,R);

% run simulation
for i = 2:nstep
    u(:,i) = -L*x(:,i-1); % calculate motor command
    x(:,i) = Ad*x(:,i-1) + Bd*u(:,i); % calculate new state
%     if i > nstep*0.5 % if desired, change target location
%         x(4,i) = 0.05;
%     end
end

figure
subplot(1,2,1)
plot(0:delt:delt*nstep-delt,x(1,:),'LineWidth',1.5)
ylabel('Position (m)')
xlabel('Time (s)')
% xlim([0 1.5])

subplot(1,2,2)
plot(0:delt:delt*nstep-delt,x(2,:),'LineWidth',1.5)
ylabel('Velocity (m/s)')
xlabel('Time (s)')
% xlim([0 1.5])
% subplot(1,3,3)
% plot(0:delt:delt*nstep-delt,u)
% hold on
% title('Motor Commands')
% xlabel('Time (s)')
% xlim([0 1])