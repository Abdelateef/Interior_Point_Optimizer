function [x0,y0,s0]=Initial_Points(A,b,c)
% intial values for the intial points 
At=transpose(A);
AAt_inv=inv(A*At);
x_hat=At*AAt_inv*b;
lamda_hat=AAt_inv*A*c;
s_hat=c-At*lamda_hat;
% eleminate the nonpositive components
sgma_x1=max((-3/2)*min(x_hat),0);
sgma_s1=max((-3/2)*min(s_hat),0);
e=ones(size(x_hat));
et=transpose(e);
x_hat1=x_hat+sgma_x1*e;
s_hat1=s_hat+sgma_s1*e;
%check that (xhat,shat)>=0
xhat1t=transpose(x_hat1);
sgma_x2=.5*((xhat1t*s_hat1)/(et*s_hat1));
sgma_s2=.5*((xhat1t*s_hat1)/(et*x_hat1));
%final values for the intial points
x0=x_hat1+sgma_x2*e;
y0=lamda_hat;
s0=s_hat1+sgma_s2*e;
end
