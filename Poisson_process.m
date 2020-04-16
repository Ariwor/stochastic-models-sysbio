%% Simulate P paths of a Poisson process

P=5; % Number of paths
lambda=5;      
T_final=10;

% Generate paths
figure();
hold on

for j=1:P
    T(1)=-log(rand)./lambda;
    i=1;
while T(i) < T_final
    T(i+1)=T(i)-log(rand)./lambda;
    i=i+1;
end
    T(i)=T_final;
    stairs(T(1:i), 0:(i-1));
end

%% Simulate P paths of a Poisson process and plot realizations of the random variable N

P=1000; % Number of paths
lambda=5;      
T_final=5;
N = zeros(1,1000);
% Generate paths
for j=1:P
    T(1)=-log(rand)./lambda;
    i=1;
while T(i) < T_final
    T(i+1)=T(i)-log(rand)./lambda;
    i=i+1;
end
    T(i)=T_final;
    N(j)=i;
end

% Plotting
figure();
tiledlayout(2,1)
nexttile
edges = [0:1:30];
h1 = histogram(N, edges, 'Normalization','probability', 'FaceColor','#77AC30');
x = [0:1:30];
pmf_poisson = poisspdf(x, 25);
nexttile
scatter(x,pmf_poisson, 'r');