clc;clear;
stocks = hist_stock_data('01012006','01012015','S&P-100.txt');

Ncomp = size(stocks,2);
Ndays = size(stocks(1).Date,1);
stocks_close = zeros(Ndays,Ncomp);

for i = 1:Ncomp
    stocks_close(:,i) = stocks(i).Close(1:Ndays);
end
save stockdata.mat
save('stocks','stocks_close')
