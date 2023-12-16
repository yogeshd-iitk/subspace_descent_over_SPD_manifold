function[data,no_comp,f_user]=user_input(W)
%% Data for algorithm
    rng("default")  % Reproducibility for Random Number Generator
    temp_C=randn(W,W);
    data.C_p=temp_C*(temp_C');
    temp_D=randn(W,W);
    data.D_q=temp_D*(temp_D');
    data.logdt=input('1 if logdet term is present 0 otheriwise=');
    temp_A=randn(W,W);
    data.A_r=temp_A*(temp_A'); 
    temp_H=randn(W,W);
    data.H_r=temp_H*(temp_H'); 
    temp_F=randn(W,W);
    data.F_s=temp_F*(temp_F'); 
    temp_G=randn(W,W);
    data.G_s=temp_G*(temp_G'); 
    %
    % temp_P=randn(W,W); % randi([1,10],W,W);
    % data.P_m=temp_P*(temp_P'); %temp_P*(temp_P')/(W^2); 
    % temp_Q=randn(W,W); % randi([1,10],W,W);
    % data.Q_m=temp_Q*(temp_Q'); % temp_Q*(temp_Q')/(W^2); 
    %
%% Number of components
[no_comp]=number_of_comp(data);
%
%%  Symbolic variables of Function
% variables
    syms a [no_comp.P 1] real
    syms b [no_comp.Q 1] real
    syms k [no_comp.logdt 1] real
    syms r [no_comp.R 1] real 
    syms s [no_comp.S 1] real 
    syms t [no_comp.M 1] real 
%
%% functional form
f_user=(sum(a)+sum(b)-sum(k)+sum(r)+sum(s)+sum(t));
%
end 