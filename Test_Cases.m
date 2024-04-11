clear
clc
%% Example 1
c=[-2;-3;0;0];
A=[2 1 1 0;1 2 0 1];
b=[4;5];
alpha=0.8;
beta=0.3;
Tol=0.01;
Central_Path(A,b,c,alpha,beta,Tol,'Fixed');
Central_Path(A,b,c,alpha,beta,Tol,'Adaptive');
Mehrotra(A,b,c,beta,Tol);
%% Example 2
% c=[-1.1;-1;0];
% A=[1 1 1];
% b=[6];
% alpha=0.5;
% beta=0.5;
% Tol=0.01;
% Central_Path(A,b,c,alpha,beta,Tol,'Fixed');
% Central_Path(A,b,c,alpha,beta,Tol,'Adaptive');
% Mehrotra(A,b,c,beta,Tol);