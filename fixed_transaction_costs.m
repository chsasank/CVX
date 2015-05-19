%% Fixed Transaction Costs
% Author: Sasank, Ashok
% Load stocks data for S&P
clc;clear;addpath('combn')
stocks_close = load('stocks.mat');
stocks_close = stocks_close.stocks_close;
n = 100;
exhaustive = false; %exponential time!
sd_sweep = 0.1:0.05:0.3;

%% Estimate means and variances of returns
t = 20;

p_now = stocks_close;
p_next = circshift(p_now,t);

p_now = p_now(t+1:end,:);
p_next = p_next(t+1:end,:);

returns = p_next./p_now;
a_bar = mean(returns)';
Sigma = cov(returns);
% a_bar = a_bar(20:20+n-1);
% Sigma = cov(returns(:,20:20+n-1));

%% Parameters

% w = 1/10*ones(n,1); % Initial portofolio holding
% alpha = 0.01*ones(n,1); % linear transacion costs
% beta = 0.001*ones(n,1); % fixed transaction costs
% s = 0.005*ones(n,1); % shorselling constraint

w = 1/100*ones(100,1); % Initial portofolio holding
alpha = 0.001*ones(100,1); % positive transacion costs
beta = 0.0005*ones(100,1); % negative transacion costs
s = 0.000*ones(100,1); % shorselling constraint
%% Solution without regard to fixed costs
w_1 = zeros(numel(sd_sweep),n);
mu_1 = zeros(numel(sd_sweep),1);
sd_1 = zeros(numel(sd_sweep),1);

delta = 1e-4;
h = waitbar(0,'Computing solution without regard to fixed costs');
for i = 1:numel(sd_sweep)
    cvx_begin quiet
        variables x_plus(n) x_minus(n)
        maximize(a_bar'*(w + x_plus - x_minus))
        subject to

        ones(n,1)'*(x_plus - x_minus) + alpha'*x_plus + alpha'*x_minus <= 0 % self financing constraint

        (w + x_plus - x_minus)'*Sigma*(w + x_plus - x_minus) <= sd_sweep(i)^2 %standard deviation constraint

        x_plus >= zeros(n,1)
        x_minus >= zeros(n,1)
        w + x_plus - x_minus >= s % shorselling constraint
    
    cvx_end
    
    w_new_i = w + x_plus- x_minus;
    w_1(i,:) = w_new_i;
    mu_1(i) = a_bar'*w_new_i;
    sd_1(i) = sqrt(w_new_i'*Sigma*w_new_i);
    
%     transacted = (abs(x_plus- x_minus) > delta);
%     extraCost = ones(n,1)'*(x_plus - x_minus) + beta'*transacted +alpha'*x_plus + alpha'*x_minus;
%     mu_1(i) = mu_1(i) - extraCost;
    waitbar(i/numel(sd_sweep),h)
end

%% Convex relaxation
w_2 = zeros(numel(sd_sweep),n);
x_plus_2 = zeros(numel(sd_sweep),n);
x_minus_2 = zeros(numel(sd_sweep),n);
mu_2 = zeros(numel(sd_sweep),1);
sd_2 = zeros(numel(sd_sweep),1);

waitbar(0,h,'Computing solution to upper bound/convexified problem');

for i = 1:numel(sd_sweep)
    
    ub = sd_sweep(i)*sqrt(diag(inv(Sigma))) - w;
    lb = -sd_sweep(i)*sqrt(diag(inv(Sigma))) - w;
    
    plus_slope = (beta./ub + alpha);
    minus_slope = (-beta./lb + alpha);
    
    cvx_begin quiet
        variables x_plus(n) x_minus(n)
        maximize(a_bar'*(w + x_plus - x_minus))
        subject to

        ones(n,1)'*(x_plus - x_minus) + plus_slope'*x_plus + minus_slope'*x_minus <= 0 % self financing constraint

        (w + x_plus - x_minus)'*Sigma*(w + x_plus - x_minus) <= sd_sweep(i)^2 %standard deviation constraint
         
        x_plus >= zeros(n,1) 
        x_plus <= ub
        
        x_minus>= zeros(n,1)
        x_minus <= -lb
        w + x_plus - x_minus >= s % shorselling constraints
    
    cvx_end
    
    x_plus_2(i,:) = x_plus;
    x_minus_2(i,:) = x_minus;
    
    w_new_i = w + x_plus- x_minus;
    w_2(i,:) = w_new_i;
    mu_2(i) = a_bar'*w_new_i;
    sd_2(i) = sqrt(w_new_i'*Sigma*w_new_i);
    waitbar(i/numel(sd_sweep),h)
end

%% Heuristic
w_3 = zeros(numel(sd_sweep),n);
mu_3 = zeros(numel(sd_sweep),1);
sd_3 = zeros(numel(sd_sweep),1);

waitbar(0,h,'Computing solution using heuristic');
delta = 1e-4;
maxIter = 20;

for i = 1:numel(sd_sweep)

    x_k = (x_plus_2(i,:) - x_minus_2(i,:))'; % initialize with solution from convexified problem

    for k = 1:maxIter   
        modified_slope = (beta./(abs(x_k)+delta) + alpha);

        cvx_begin quiet
            variables x_plus(n) x_minus(n)
            maximize(a_bar'*(w + x_plus - x_minus))
            subject to

            ones(n,1)'*(x_plus - x_minus) + modified_slope'*x_plus + modified_slope'*x_minus <= 0 % self financing constraint

            (w + x_plus - x_minus)'*Sigma*(w + x_plus - x_minus) <= sd_sweep(i)^2 %standard deviation constraint

            x_plus >= zeros(n,1) 
            x_minus >= zeros(n,1)  
            w + x_plus - x_minus >= s % shorselling constraints

        cvx_end
        
        %check for convergence
        if norm(x_plus - x_minus - x_k)/norm(x_k) < 1e-6
            break;
        end
        
        x_k = x_plus - x_minus;
    end
    w_new_i = w + x_plus- x_minus;
    w_3(i,:) = w_new_i;
    mu_3(i) = a_bar'*w_new_i;
    sd_3(i) = sqrt(w_new_i'*Sigma*w_new_i);
    waitbar(i/numel(sd_sweep),h)
end
%% Exhaustive search
if(exhaustive)
tic
w_4 = zeros(numel(sd_sweep),n);
mu_4 = zeros(numel(sd_sweep),1);
sd_4 = zeros(numel(sd_sweep),1);

waitbar(0,h,'Computing using Exhaustive search');
M = combn([0 1],n); % generate all possible 0,1 vectors

for i = 1:numel(sd_sweep)
    w_best = zeros(n,1);
    mu_best = 0;
    for j = 1:10%size(M,1)       
        %M(j,i) specifies which asset has to be zero
        
        cvx_begin quiet
            variables x_plus(n) x_minus(n)
            maximize(a_bar'*(w + x_plus - x_minus))
            subject to

            ones(n,1)'*(x_plus - x_minus) + (~M(j,:))*beta +alpha'*x_plus + alpha'*x_minus <= 0 % self financing constraint

            (w + x_plus - x_minus)'*Sigma*(w + x_plus - x_minus) <= sd_sweep(i)^2 %standard deviation constraint

            x_plus >= zeros(n,1)
            x_minus >= zeros(n,1)
            w + x_plus - x_minus >= s % shorselling constraint
            
            M(j,:)*x_plus == 0
            M(j,:)*x_minus == 0

        cvx_end
        
        w_current = w + x_plus- x_minus;
        mu_current = a_bar'*w_current;
        if mu_current >= mu_best
            w_best = w_current;
            mu_best = mu_current;
        end
    end
    w_4(i,:) = w_best;
    mu_4(i) = a_bar'*w_best;
    sd_4(i) = sqrt(w_best'*Sigma*w_best);
    waitbar(i/numel(sd_sweep),h)
end
toc
end
close(h)
%% Results
figure
plot(sd_1,mu_1/sum(w),'c')
hold on
plot(sd_2,mu_2/sum(w),'r')
plot(sd_3,mu_3/sum(w),'g')
if(exhaustive)
    plot(sd_4,mu_4/sum(w),'b')
    legend('without regard to fixed costs','upper bound/convexifaction','heuristic','exhaustive')
else
    legend('without regard to fixed costs','upper bound/convexifaction','heuristic')
end

%%
w_new = w+x_plus- x_minus;
d= [w_new , x_plus-x_minus, a_bar-1, diag(Sigma)]*100;
d = d./repmat(sqrt(sum(d.^2)),n,1); %normalize so that can be plotted in one plot

% mu_new = a_bar'*w_new;
% sd_new = sqrt(w_new'*Sigma*w_new);
% z = 0:0.001:2;
% figure
% plot(z,normcdf(z,mu_new,sd_new))

figure
plot(d(:,2))
hold on
plot(d(:,3),'r')
plot(d(:,4),'g:')
legend('Transactions x','mean returns','variance in return')
xlabel('Stock')