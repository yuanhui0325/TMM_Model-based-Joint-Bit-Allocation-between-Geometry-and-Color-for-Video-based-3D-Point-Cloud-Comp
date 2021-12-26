function [c,ceq] = fmincon_rate(x,n)
c(1) = n(1).*power(x(1),n(2))+n(3).*power(x(2),n(4))-n(5);
ceq = [];