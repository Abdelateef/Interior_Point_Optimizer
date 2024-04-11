function Central_Path(A,b,c,alpha,beta,Tol,alpha_type)
[x,y,s]=Initial_Points(A,b,c); %initals 
X_points=[x];
S_points=[s];
[m,n]=size(A);
%% J matrix componants
J11=zeros(n);
At=transpose(A);
I= eye(n);
J22=zeros(m,n);
J23=zeros(m);
S=diag(s);
J32=zeros(n,m);
X=diag(x);
%% f componants
e=ones(n,1);
Tao=transpose(x)*s;
Tao_array=[Tao];
mu=(Tao/n);
mu_array=[mu];
rc=At*y+s-c;
rb=A*x-b;
rxs=x.*s- mu.*beta.*e;
k =1;
No_iteration=[k];
while((mu)>=Tol)
I= eye(n);
J=[J11, At, I; 
    A, J22, J23;
    S, J32, X];
f=[-rc;-rb;-rxs];
%calculate J by normal inverse
deltas=inv(J)*f;
deltas_x=deltas(1:n);
[nn,mm]=size(y);
deltas_y=deltas(n+1:nn+n);
deltas_s=deltas(nn+n+1:nn+2*n);

%% update
if strcmp(alpha_type, 'Fixed')
    prim_alpha=alpha;
    dual_alpha=alpha;
elseif strcmp(alpha_type, 'Adaptive')
    x_ratio = min(x./abs(deltas_x));
    s_ratio=min(s./abs(deltas_s));
    prim_alpha=min(1,x_ratio);
    dual_alpha=min(1,s_ratio);
else 
    disp('you should enter the alpha_type Fixed or Adaptive');
end 
x=x+prim_alpha*deltas_x;
s=s+dual_alpha*deltas_s;
if (all(x>=0) && all(s>=0)) %check the entring values of s and x are postive
    y=y+dual_alpha*deltas_y;
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
    rxs=x.*s- mu.*beta.*e;
    k =k+1;
    No_iteration=[No_iteration,k];
else 
    break;
end

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
title([alpha_type,' Step Size','- Central Path Method','- Feasible Region']);
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
title([alpha_type,' Step Size','- Central Path Method',' - Complementarity Condition']);

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
title([alpha_type,' Step Size','- Central Path Method',' - Objective Function']);
end