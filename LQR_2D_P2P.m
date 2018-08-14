% Simulation of a point to point reach with different degrees of rotation.
% This script produces 2D position plots of the three angles of rotation
% the user selects. Parameters to adjust before starting: delt, nstep,
% distance, ang, Q, R.

clear all;
rand('seed',1);

% adjust time step 
delt = 0.00001; % time step length in secs
nstep = ceil(0.7/delt);%5000; % number of time steps
% delt = 0.01;
% nstep = 100;
distance = -8; % horizontal distance of hand from zero
ang = [0 20 30 45 75 90]; % six different angles of rotation to simulate

% values for Q and R taken from Qian infinite horizon model
Q = diag([1 0.01 0 1 0.01 0])*0.001; % accuracy cost- default is [1 0.01 0]
R = eye(2)*0.0001; % effort cost- default is 0.0001

z = 6; % number of simulations
order = 6; % order of the system

% parameters for A and B matrices
t1 = 0.224;
t2 = 0.013;
t3 = 0.004;
k = 0; %1 to include spring
b = t1 + t2;
m = t1*t2;
r = t3;

% generate A and B matrices in discrete time formulation
y = [0 0 0];
A = [0 1 0 y; -k/m -b/m 1/m y; 0 0 -1/r y; y 0 1 0; y -k/m -b/m 1/m; y 0 0 -1/r];
B = [0 0 1/r y; y 0 0 1/r]';
% Ad = eye(order) + delt*A;
Ad = expm(delt*A);
Bd = delt*B;

[P,eig,L] = care(A,B,Q,R);

xt = zeros(order,nstep,z);
xt(1,1,:) = distance; % initialize state variables; x position
xt(4,1,:) = 0; % y position -> keep at 0
u = zeros(2,nstep,z); % movement commands
L2 = L;

for j = 1:z
    % calculate rotation matrix and apply to L
    R2 = rotz(ang(j));
    R2 = R2(1:2,1:2);
    L = R2*L2;
    
    for i = 2:nstep
        u(:,i,j) = -L*xt(:,i-1,j);
        xt(:,i,j) = Ad*xt(:,i-1,j) + Bd*u(:,i,j);
    end
end

% plot 2D hand trajectories for each rotation
a = [distance 2 -1 4];
figure;
subplot(2,3,1)
plot(xt(1,:,1),xt(4,:,1),'LineWidth',3);
axis(a);
title([num2str(ang(1)),' degrees'],'FontSize',15)
ylabel('Position','FontSize',15);

subplot(2,3,2)
plot(xt(1,:,2),xt(4,:,2),'LineWidth',3);
axis(a);
title([num2str(ang(2)),' degrees'],'FontSize',15);

subplot(2,3,3)
plot(xt(1,:,3),xt(4,:,3),'LineWidth',3);
axis(a);
title([num2str(ang(3)),' degrees'],'FontSize',15);

subplot(2,3,4)
plot(xt(1,:,4),xt(4,:,4),'LineWidth',3);
axis(a);
title([num2str(ang(4)),' degrees'],'FontSize',15)
ylabel('Position','FontSize',15);

subplot(2,3,5)
plot(xt(1,:,5),xt(4,:,5),'LineWidth',3);
axis(a);
title([num2str(ang(5)),' degrees'],'FontSize',15);
xlabel('Position','FontSize',15);

subplot(2,3,6)
plot(xt(1,:,6),xt(4,:,6),'LineWidth',3);
axis(a);
title([num2str(ang(6)),' degrees'],'FontSize',15);