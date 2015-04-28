%% Load stocks data for S&P
clc;clear;
stocks_close = load('stocks.mat');
stocks_close = stocks_close.stocks_close;

%% Find 20 day returns
n = 20;

p_now = stocks_close;
p_next = circshift(p_now,n);

p_now = p_now(n+1:end,:);
p_next = p_next(n+1:end,:);

returns = p_next./p_now;
%% Estimate means and variances of returns
a = mean(returns);
sigma = cov(returns);