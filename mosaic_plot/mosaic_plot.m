% Make a mosaic plot (please also see
% http://en.wikipedia.org/wiki/Mosaic_plot)
% N.B. data should be non-negative
%
% outputs xm,ym are x and y components of centers of boxes. Ther are useful
% for labeling the plot. 
%
% E.g., making a mosaic plot for random data and mark the percentage of box as
% label
% 
% data=rand(3,4);
% [xm,ym]=mosaic_plot(data);
% multi_text(xm(:),ym(:),form_percentage_strings_from_array(data(:)));
 

function [xm,ym]=mosaic_plot(data)

if min(data(:)) < 0
    error('data has to be non-negative');
end

xs=sum(data);
xs=xs/sum(xs);
ys=data*diag(1./sum(data));

gap=min([xs(:); ys(:)])/4;
gap=min([gap,0.01]);

xs=[0 cumsum(xs)];
ys=cumsum(ys);
ys=[zeros(1,size(ys,2)); ys];

for id=1:length(xs)-1
   plot_rectangles(xs(id),xs(id+1),ys(:,id),gap); 
end

ym=diff(ys)/2+ys(1:end-1,:);
xm=diff(xs)/2+xs(1:end-1);
xm=ones(size(ym,1),1)*xm;

