function nice_plot(figNum, myXLABEL, myYLABEL, square)

mat = get(0,'ScreenSize');
w = mat(3)*0.7; l = mat(4)*0.7;

if square == false
    width = 800;
    height = 400;
else
    width = 600;
    height = 600;
end


figure(figNum);
set(gca,'FontSize', 20);

h = xlabel(myXLABEL);
set(h, 'Interpreter','Latex', 'FontSize', 25);

h = ylabel(myYLABEL);
set(h, 'Interpreter','Latex', 'FontSize', 25);

box off;
set(gcf,'color','w');

p = get(gcf,'Position');
set(gcf,'Position',[p(1),p(2)+500,width,height]);


end