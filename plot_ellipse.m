function plot_ellipse(mu, Sigma)
% plots an ellipse around mu, illustrating covariance matrix Sigma
N=100;
t = linspace(0,2*pi,N);
pp = [cos(t); sin(t)];

[V D] = eig(Sigma);
VV = V*sqrt(D);

p = VV*pp + repmat(mu,1,N);
plot(p(1,:),p(2,:),'b');