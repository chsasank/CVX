function multi_text(xs,ys,texts) % put text several times

for id=1:length(xs(:))
    h=text(xs(id),ys(id),texts{id});
    set(h,'horizontalalignment','center');
end