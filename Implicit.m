%Implicit Scheme to solve Black Scholes PDE with up and out call option
clear all, close all
global T alpha
%% Parameters
K = 90; % Strike Price
So=100; %Spot Price 
r = 0.03;% Interest rate
q = 0.05;% Dividend yield
B = 130;% Barrier Level 
alpha=0.35;% Exponent in local volativity function
T = 0.5;% Time to maturity
M = 10000; %Number of Time Steps
N = 10; % Number of Division

% Set the minimal and maximal stock prices
Smin = 0;
Smax = B;

%% Numerical Discretization setting
% Setup our grid in stock price direction
S1 = linspace(Smin,Smax,N+1)';
dS = S1(2) - S1(1); % Grid cell size
S = S1(2:N); % S stores all the prices except boundary points

% Setup our grid in time direction
tau = linspace(0,T,M+2); % time values evaluated
dtau = tau(2) - tau(1); % Time Step magnitude

%Solution matrix to store value of V(t,S) at all grid points
V = zeros(N-1,M+1);
% Vector to store the option prices at time k
% For initilization the up - and-out call option was used payoff
% V= max(S-K,0), if S<B
%    0,         otherwise
%
V(:,1)= max([S-K,zeros(size(S))],[],2);


%% Solve BS PDE (IMPLICIT)
for k=1:M+1
    % lower, main and upper diagonal of the tridiagonal matrix for the implicit
    % finite difference scheme
    l=-(0.5.*(dtau/(dS^2))*(sigma(tau(k),S(:)).^2).*(S(:).^2)-0.5*(dtau/dS)*(r-q).*S(:));
    d=-(-(dtau/(dS^2))*(sigma(tau(k),S(:)).^2).*(S(:).^2)-r*dtau-1);
    u=-(0.5*(dtau/(dS^2))*(sigma(tau(k),S(:)).^2).*(S(:).^2) +0.5*(dtau/dS)*(r-q).*S(:));
    
    V(:,k+1)=tridiag (d,u,l,V(:,k));
end

%% Results
% Interpolation to find the put price when S0 = 100.0
call_fdm=interp1(S,V(:,end),So);

% % Plot of FDM BS function
figure()
plot(S,V(:,end),'LineWidth',2)
title('Implicit - BS formula')
xlabel('Stock price')
ylabel('Call price')
legend('Implicit','Location','SouthEast')

% % 3D surface Plot of payoff price FDM BS function
figure()
surf(tau,S,V,'edgecolor','none')
title('Implicit - BS formula')
xlabel('Time (years)')
ylabel('Stock price')
zlabel('Put price')
legend('Implicit','Location','SouthEast')

%% Functions
%local volatibity function
function resp=sigma(ti,Sn)
    global T alpha
    resp=0.25.*exp(T-ti).*(100./Sn).^alpha;
end

%Function to solve tridiagonal matrix linear system Ax=b
function x = tridiag (D,U,L,B)
  n = length(B);
  x = zeros(n,1);

  %make A upper diagonal
  for j=2:n
      D(j) = D(j) - L(j)*U(j-1) / D(j-1);
      B(j) = B(j) - L(j)*B(j-1) / D(j-1);
  end
  %Back substitution
  x(n) = B(n) / D(n);
  for j=n-1:-1:1
      x(j) = (B(j) - U(j)*x(j+1)) / D(j);
  end
end