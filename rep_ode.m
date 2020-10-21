function xdot = rep_ode(t, x, A, Q, lambda)

    f1 = ((A(1,1)*x(1)) + (A(1,2)*x(2)) + (A(1,3)*x(3)) + (A(1,4)*x(4)))/sum(x(1:4));
    f2 = ((A(2,1)*x(1)) + (A(2,2)*x(2)) + (A(2,3)*x(3)) + (A(2,4)*x(4)))/sum(x(1:4));
    f3 = ((A(3,1)*x(1)) + (A(3,2)*x(2)) + (A(3,3)*x(3)) + (A(3,4)*x(4)))/sum(x(1:4));
    f4 = ((A(4,1)*x(1)) + (A(4,2)*x(2)) + (A(4,3)*x(3)) + (A(4,4)*x(4)))/sum(x(1:4));
    phi0 = x(1)*f1 + x(2)*f2 + x(3)*f3 + x(4)*f4; % dot product
    
    % yes mutations
    xdot(1,1) = x(1)*f1*Q(1,1) + x(2)*f2*Q(2,1) + x(3)*f3*Q(3,1) + x(4)*f4*Q(4,1) - phi0*x(1);
    xdot(2,1) = x(1)*f1*Q(1,2) + x(2)*f2*Q(2,2) + x(3)*f3*Q(3,2) + x(4)*f4*Q(4,2) - phi0*x(2);
    xdot(3,1) = x(1)*f1*Q(1,3) + x(2)*f2*Q(2,3) + x(3)*f3*Q(3,3) + x(4)*f4*Q(4,3) - phi0*x(3);
    xdot(4,1) = x(1)*f1*Q(1,4) + x(2)*f2*Q(2,4) + x(3)*f3*Q(3,4) + x(4)*f4*Q(4,4) - phi0*x(4);
    
    vvv = x(5) + x(6) + x(7) + x(8);
    
    xdot(5,1) = (lambda + f1)*x(1)*vvv;
    xdot(6,1) = (lambda + f2)*x(2)*vvv;
    xdot(7,1) = (lambda + f3)*x(3)*vvv;
    xdot(8,1) = (lambda + f4)*x(4)*vvv;
    
end