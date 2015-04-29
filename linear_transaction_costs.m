%% Linear Transaction Costs
% Load stocks data for S&P
clc;clear;
stocks_close = load('stocks.mat');
stocks_close = stocks_close.stocks_close;
n = size(stocks_close,2);
%% Find 20 day returns
t = 20;

p_now = stocks_close;
p_next = circshift(p_now,t);

p_now = p_now(t+1:end,:);
p_next = p_next(t+1:end,:);

returns = p_next./p_now;
%% Estimate means and variances of returns
a_bar = mean(returns)';
sigma = cov(returns);

%% Parameters
w = 1/100*ones(100,1); % Initial portofolio holding
alpha_plus = 0.01*ones(100,1); % positive transacion costs
alpha_neg = 0.01*ones(100,1); % negative transacion costs
s = 0.005*ones(100,1); % shorselling constraint

eta_1 = 0.80; % confidence for shortfall 1
W_low_1 = 0.9; % bad return

eta_2 = 0.97; % confidence
W_low_2 = 0.7; % disastrous return

%% Optimization Problem
cvx_begin %quiet
    variables x_plus(n) x_neg(n)
    maximize(a_bar'*(w + x_plus - x_neg))
    subject to 
        ones(n,1)'*(x_plus - x_neg) + alpha_plus'*x_plus + alpha_neg'*x_neg <= 0 % self financing constraint
        x_plus >= zeros(n,1)
        x_neg >= zeros(n,1)
        w + x_plus - x_neg >= s % shorselling constraint
        
        norminv(eta_1)*norm(sigma^0.5*(w + x_plus - x_neg)) <= a_bar'*(w + x_plus - x_neg) - W_low_1 % shortfall constraint 1
        norminv(eta_2)*norm(sigma^0.5*(w + x_plus - x_neg)) <= a_bar'*(w + x_plus - x_neg) - W_low_2 % shortfall constraint 2

cvx_end

%% Results
w_new = w+x_plus- x_neg;
d= [w_new , x_plus-x_neg, a_bar-1, diag(sigma)]*100

mu_new = a_bar'*w_new;
sd_new = sqrt(w_new'*sigma*w_new);
z = 0:0.001:2;
plot(z,normcdf(z,mu_new,sd_new))

figure
plot(d(:,2))
hold on
plot(d(:,3),'r')
plot(d(:,4),'g:')
legend('Transactions x','mean returns','variance in return')
xlabel('Stock')