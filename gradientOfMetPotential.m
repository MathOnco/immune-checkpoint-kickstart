function c = gradientOfMetPotential(b3,b4)
    b1 = 0; b2 = 0;
    xlim = 101;
    ylim = 101;

    c = zeros(xlim,ylim);

    bMax = b3 + b4;
    for i = 1:1:(xlim)
        for j = 1:1:(ylim)

            xp = (i-1)/100;
            yp = (j-1)/100;

            x0 = b1 + (b3 - b1)*xp;
            y0 = b2 + (b4 - b2)*yp;
            c(i,j) = (x0+y0)/bMax;
        end
    end
end