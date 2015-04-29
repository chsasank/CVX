%% Linear Transaction Costs
% Load stocks data for S&P
clc;clear;
stocks_close = load('stocks.mat');
stocks_close = stocks_close.stocks_close;
n = 10;
%% Find 20 day returns
t = 20;

p_now = stocks_close;
p_next = circshift(p_now,t);

p_now = p_now(t+1:end,:);
p_next = p_next(t+1:end,:);

returns = p_next./p_now;
%% Estimate means and variances of returns
a_bar = mean(returns)';
Sigma = cov(returns);

[~,sortIndex] = sort(a_bar,'descend');
a_bar = a_bar(sortIndex(1:n));
Sigma = cov(returns(:,sortIndex(1:n)));
%% Parameters

w = 1/10*ones(10,1); % Initial portofolio holding
alpha = 0.01*ones(10,1); % linear transacion costs
beta = 0.01*ones(10,1); % fixed transaction costs
s = 0.05*ones(10,1); % shorselling constraint

sd_sweep = 0.1:0.005:0.3;

%% Solution without regard to fixed costs
w_1 = zeros(numel(sd_sweep),n);
mu_1 = zeros(numel(sd_sweep),1);
sd_1 = zeros(numel(sd_sweep),1);

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
    waitbar(i/numel(sd_sweep),h)
end

%% Convex relaxation
w_2 = zeros(numel(sd_sweep),n);
x_plus_2 = zeros(numel(sd_sweep),n);
x_minus_2 = zeros(numel(sd_sweep),n);
mu_2 = zeros(numel(sd_sweep),1);
sd_2 = zeros(numel(sd_sweep),1);

waitbar(0,h,'Computing solution to upper bound/convexifaction');

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
         
        zeros(n,1) <= x_plus <= ub
        zeros(n,1) <= x_minus <= -lb
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

waitbar(0,h,'Computing solution for heuristic');
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

            zeros(n,1) <= x_plus <= ub
            zeros(n,1) <= x_minus <= -lb
            w + x_plus - x_minus >= s % shorselling constraints

        cvx_end
        
       
        if norm(x_plus - x_minus - x_k)/norm(x_k) < 1e-4
            break;
        end
        
        x_k = x_plus - x_minus;
    end
    k
    w_new_i = w + x_plus- x_minus;
    w_3(i,:) = w_new_i;
    mu_3(i) = a_bar'*w_new_i;
    sd_3(i) = sqrt(w_new_i'*Sigma*w_new_i);
    waitbar(i/numel(sd_sweep),h)
end
%% Results
figure
plot(sd_1,mu_1)
hold on
plot(sd_2,mu_2,'r')
plot(sd_3,mu_3,'g')
legend('without regard to fixed costs','upper bound/convexifaction','heuristic')
