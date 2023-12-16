function[data,no_comp,f_user]=user_input(W)
%% Data for algorithm
    rng("default")  % Reproducibility for Random Number Generator
    P=4;
    for p=1:P
    temp_C=randn(W,W);
    data.C_p(:,:,p)=temp_C*(temp_C');
    end
    Q=4;
    for q=1:Q
    temp_D=randn(W,W);
    data.D_q(:,:,q)=temp_D*(temp_D');
    end
    data.logdt=input('1 if logdet term is present 0 otheriwise=');
    R=5;
    for r=1:R
    temp_A=randn(W,W);
    data.A_r(:,:,r)=temp_A*(temp_A'); 
    temp_H=randn(W,W);
    data.H_r(:,:,r)=temp_H*(temp_H'); 
    end
    S=6;
    for s=1:S
        temp_F=randn(W,W);
        data.F_s(:,:,s)=temp_F*(temp_F'); 
        temp_G=randn(W,W);
        data.G_s(:,:,s)=temp_G*(temp_G'); 
    end
    %
    % M=8;
    % for t=1:M
    %     temp_P=randn(W,W); 
    %     data.P_m(:,:,t)=temp_P*(temp_P'); 
    %     temp_Q=randn(W,W); 
    %     data.Q_m(:,:,t)=temp_Q*(temp_Q'); 
    % end
    %
%% Number of components
[no_comp]=number_of_comp(data);
%
%%  Symbolic variables of Function
% variables
    syms p [no_comp.P 1] real
    syms q [no_comp.Q 1] real
    syms k [no_comp.logdt 1] real
    syms r [no_comp.R 1] real 
    syms s [no_comp.S 1] real 
    syms t [no_comp.M 1] real 
%
%% functional form
f_user=(sum(p)+sum(q)-sum(k)+sum(r)+sum(s)+sum(t));
%
end 
