function PlotMarkers(x,y)
cla;
cols = get(gca,'ColorOrder');
cols(4:6,:) = cols(1:3,:);
for i = 1:size(x,1)
    if y(i) <= 7
        col = cols(y(i),:);
    else
        col = [0 0 0];
    end
    
    plot3(x(i,1),x(i,2),x(i,3),'.','MarkerSize',20,'Color',col)
    hold on
end
view([0 1 0]);
axis equal;
end