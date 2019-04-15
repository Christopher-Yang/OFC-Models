% simulate multiple times with noise to get a noisy phaser plot
clear all
load('av_data.mat');

Nreps = 100;

L = [80 9 5];

for i=1:Nreps
    sim(i) = freq_sim_noisy(L,i);
end

%% plot trajectory of hand and target
figure(1); clf; hold on
for i=1:Nreps
    plot(sim(i).t,sim(i).target,'y')
    plot(sim(i).t,sim(i).hand,'b','LineWidth',1);
end
%title('Simulation trajectory');
leg2 = {'Target','Hand'};
legend(leg2,'FontSize',15);
xlabel('Time','FontSize',15);
ylabel('Position','FontSize',15);
xlim([0 10]);
% legend('Target','.0001x effort','.01x','Default','100x','10000x');
% legend('Target','10 ms delay','100 ms','200 ms (default)','400 ms','600 ms');

%% plot phasors
figure(4); clf; hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
Nfreq = length(sim(1).ratio);

% plot the data
%data.freqs_x = data.rot.subj17.no_rot1.freqs_x;
%plot(data.rot.avg.x_x.fft(1,:),'ro','LineWidth',1.5);
plot(data.fft(1,:),'ro','LineWidth',1.5);

for i=1:Nreps
    for j=1:1:Nfreq
        ratio(i,j) = sim(i).ratio(j);
    end
end

for i=1:length(data.freqs_x)
    f_ind(i) = find(abs(sim(1).freq - data.freqs_x(i))<.001); % find corresponding frequency for data in the simulation
end

%colmap = [linspace(
ji = [];
for j=f_ind
    plot(ratio(:,j),'b.'); hold on;
    ji = [ji j];
end
axis equal



% estimate Gaussian distribution for real and complex parts
for i=1:length(ji)
    z(1,:,i) = real(ratio(:,ji(i)));
    z(2,:,i) = imag(ratio(:,ji(i)));
    
    % get mean
    z_mean(:,i) = mean(z(:,:,i)')';
    %plot(z_mean(1,i),z_mean(2,i),'ro')
    
    % get covariance matrix
    zN(:,:,i) = z(:,:,i)-repmat(z_mean(:,i),1,Nreps);
    Sigma(:,:,i) = zN(:,:,i)*zN(:,:,i)'/Nreps;
    plot_ellipse(z_mean(:,i),Sigma(:,:,i));
    
end



%polarplot(data.rot.avg.x_x.fft(1,:),'-o','LineWidth',1.5);
% polarplot(data.rot.avg.x_x.fft(2,:),'-o','LineWidth',1.5);
% polarplot(data.rot.avg.x_x.fft(5,:),'-o','LineWidth',1.5);

% legend('Inf Horizon','X -> X data','Y -> Y data');
% legend('.0001x effort','.01x','Default','100x','10000x','X -> X data','Y -> Y data');
% legend('10 ms delay','100 ms','200 ms (default)','400 ms','600 ms','X -> X data','Y -> Y data');