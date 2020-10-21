function [score, met_over_time] = calcMetScore(ode_time,x_out,y_out,p1,p2)

    bMax = p1 + p2;
    del_t = ode_time(2) - ode_time(1);
    total_del_t = del_t*length(ode_time); % ~ 1 month.
    
    score = 0;
    met_over_time = [];
    for j = 1:1:length(ode_time)
        
        x0 = p2*x_out(j);
        y0 = p1*y_out(j);

        score = score + ((x0+y0)/bMax)*del_t/total_del_t;
        met_over_time(1,j) = score;
    end
end