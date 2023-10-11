%Monte Carlo method to solve Black Scholes PDE with up and out call option
clear all
global T alpha 
%% Parameters
K = 90; % Strike Price
So=100; %Spot Price 
r = 0.03;% Interest rate
q = 0.05;% Dividend yield
B = 130;% Barrier Level 
alpha=0.35;% Exponent in local volativity function
T = 0.5;% Time to maturity
M=1000; %Number of Time Steps
N = 150; % Number of Division
NSim = 100000;%Number of MC simulations

% Set the minimal and maximal stock prices
Smin = 0;
Smax = B;

%% Numerical Discretization setting
dtau = T/M;% Time Step magnitude

% Setup our grid in stock price direction
S1 = linspace(Smin,Smax,N+1)';
dS = S1(2) - S1(1); % Grid cell size
S = S1(2:N); % S stores all the prices except boundary points

%Solution matrix to store value of V(tau=T,S) at all grid points
V = zeros(N-1,1);

%% Solve BS MC
for i=1:N-1
S_paths = zeros(NSim,M+1);
lnS1 = zeros(NSim,M+1);% init the logspot price path
lnS1(:,1)=log(S(i)*exp(-q*T));% adjust due to dividend yield
dW = randn(NSim,M);% precompute all randoms

for k=2:M+1
    lnS1(:,k) = lnS1(:,k-1) + (r-q-0.5.*sigma(k*dtau,exp(lnS1(:,k-1))).^2)*dtau + sigma(k*dtau,exp(lnS1(:,k-1))).* sqrt(dtau).*dW(:,k-1);
end

S_paths(:,:) = exp(lnS1);
V(i)=mean(prod((S_paths(:,:)<B),2).*max(S_paths(:,end)-K,0));
end 
%% Results
% Interpolation to find the put price when S0 = 100.0
call_fdm=interp1(S,V(:),So);

% % Plot of FDM BS function
figure()
plot(S,V(:),'LineWidth',2)
title('Monte Carlo - BS formula')
xlabel('Stock price')
ylabel('Call price')
legend('MC','Location','SouthEast')

%% Functions
function resp=sigma(ti,Sn)
    global T alpha
    resp=0.25.*exp(T-ti).*(100./Sn).^alpha;
end