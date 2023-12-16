clc;
clear variables;
close all;
format long
filename = ['user_data.mat'];
load(filename,'W');
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon=[];
L_0 = eye(W);
tic
options = optimoptions("fmincon",'MaxFunctionEvaluations',6000);
L = fmincon(@user_fun_con,L_0,A,b,Aeq,beq,lb,ub,nonlcon,options);
toc
f=user_fun_con(L)