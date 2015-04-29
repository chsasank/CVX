function plot_rectangles(x1,x2,ys,gap,labels,colors)

if ~exist('colors','var')
    colors={'r','b','g','k','c','m','y'};
end

if ~exist('gap','var')
    gap=0.01;
end

g2=gap/2;
xtmp=[x1+g2, x2-g2, x2-g2, x1+g2, x1+g2];
hold on;
for id=1:length(ys)-1
   ytmp=[ys(id)+g2, ys(id)+g2, ys(id+1)-g2, ys(id+1)-g2, ys(id)+g2];
   plot(xtmp,ytmp,colors{mod(id,length(colors))+1});
   
   if exist('labels','var')
      if length(labels)>=id
         h=text((x1+x2)/2,(ys(id)+ys(id+1))/2,labels{id});
         set(h,'horizontalalignment','center');
         set(h,'verticalalignment','middle');
      end
   end
end
hold off;
    

