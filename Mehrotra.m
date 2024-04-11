function Mehrotra(A,b,c,beta,Tol)
[x,y,s]=Initial_Points(A,b,c); %initals 
X_points=[x];
S_points=[s];
[m,n]=size(A);
%% J matrix componants
At=transpose(A);
S=diag(s);
X=diag(x);
%% f componants
e=ones(n,1);
Tao=transpose(x)*s;
Tao_array=[Tao];
mu=(Tao/n);
mu_array=[mu];
rc=At*y+s-c;
rb=A*x-b;
rxs_aff=x.*s;
k =1;
No_iteration=[k];
while((mu)>=Tol)
S_inv=inv(S);
D2=S_inv*X;
AD2AT_inv=inv(A*D2*At);
%%cal aff
deltas_y_aff =AD2AT_inv*(-rb-A*S_inv*X*rc+A*S_inv*rxs_aff);
deltas_s_aff=-rc-At*deltas_y_aff;
deltas_x_aff=-S_inv*rxs_aff-X*S_inv*deltas_s_aff;

rxs=x.*s + deltas_x_aff.*deltas_s_aff.*e - mu.*beta.*e;
deltas_y =AD2AT_inv*(-rb-A*S_inv*X*rc+A*S_inv*rxs);
deltas_s=-rc-At*deltas_y;
deltas_x=-S_inv*rxs-X*S_inv*deltas_s;

%cal alpha_prime alpha_daul mu_aff
x_ratio_aff= min(x./abs(deltas_x_aff));
s_ratio_aff=min(s./abs(deltas_s_aff));
prim_alpha_aff=min(1,x_ratio_aff);
dual_alpha_aff=min(1,s_ratio_aff);
mu_aff= (transpose(x+prim_alpha_aff*deltas_x_aff)*(s+dual_alpha_aff*deltas_s_aff))/n;
beta=(mu_aff/mu)^3;

% update
x=x+prim_alpha_aff*deltas_x;
s=s+dual_alpha_aff*deltas_s;
y=y+dual_alpha_aff*deltas_y;
S_points=[S_points,s];
X_points=[X_points,x];
S=diag(s);
X=diag(x);
Tao=transpose(x)*s;
Tao_array=[Tao_array,Tao];
mu=(Tao)/n;
mu_array=[mu_array,mu];
rc=At*y+s-c;
rb=A*x-b;
rxs_aff=x.*s;
k =k+1;
No_iteration=[No_iteration,k];
end
%% ploting 
%plot fesable regoin 
x_coordinates = X_points(1, :);
y_coordinates = X_points(2, :);

% Create a figure 1  
figure;
hold on
plot(x_coordinates, y_coordinates, 'x-', 'LineWidth', 2);
% % Plot the linear constraints
for i = 1:m
    st_array = [A(i, 1), A(i, 2), b(i)];
    x2 = @(x) (st_array(3) - st_array(1) * x) / st_array(2);
    fplot(x2,[0,7],'LineWidth', 2, 'Color', 'black');
end
ylim([0,7]);
xlabel('X1');
ylabel('X2');
title(['Mehrotra Method','- Feasible Region']);
hold off;
grid on;
% Create a figure 2
S1= S_points(1, :);
S2= S_points(2, :);
X1S1=[];
X2S2=[];
[nn,mm]=size(x_coordinates);
for i=1:mm
    X1S1=[X1S1,x_coordinates(1,i)*S1(1,i)];
    X2S2=[X2S2,y_coordinates(1,i)*S2(1,i)];
end 
figure;
plot(X1S1, X2S2, 'x-', 'LineWidth', 2);
xlabel('X1S1');
ylabel('X2S2');
title(['Mehrotra Method',' - Complementarity Condition']);

%Create figure 3 
f=@(x_1,x_2) c(1)*x_1+c(2)*x_2;
f_values=[];
for i=1:mm
    f_values=[f_values,-1*f(x_coordinates(1,i),y_coordinates(1,i))];
end 
figure;
plot(No_iteration, f_values, 'x-', 'LineWidth', 2, 'Color', 'red');
xlabel('iteration');
ylabel('f_max');
title([ 'Mehrotra methods',' - Objective Function']);
end