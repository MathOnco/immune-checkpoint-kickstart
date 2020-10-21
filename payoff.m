function A = payoff(kT,d1,d2,b2,c1,c2)
    A = [(-kT - d1), (-kT - d1 + b2 ), (-kT - d1), (-kT - d1 + b2 );
         (-d1 - d2*kT - c2), (-d1 - d2*kT - c2), (-d1 -d2*kT - c2), (-d1 -d2*kT - c2);
         (-kT - c1), (-kT - c1 + b2), (-kT - c1), (-kT - c1 + b2 );
         (-d2*kT - c1 - c2), (-c1 - d2*kT - c2), (-c2 - c1 -d2*kT), (-c2 - c1 - d2*kT)];
end

