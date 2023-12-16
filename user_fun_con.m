function[f]=user_fun_con(L)
% X=eye(W);
filename = ['user_data.mat'];
load(filename,'W','data','no_comp');
X=L*L'+eps*eye(W);
%% g function block calculations   
            if no_comp.P~=0
                fun_block=arrayfun(@(i) trace(data.C_p(:,:,i)/X), 1:no_comp.P); 
            else
                fun_block=[];
            end
            if no_comp.Q~=0
                fun_block=[fun_block (arrayfun(@(i)  trace(data.D_q(:,:,i)*X), 1:no_comp.Q))];
            end
            if no_comp.logdt~=0
                fun_block= [fun_block log(det(X))]; 
            end
            if no_comp.R~=0
                fun_block= [fun_block (arrayfun(@(i) trace(data.A_r(:,:,i)*X*data.H_r(:,:,i)*X), 1:no_comp.R))]; 
            end
            if no_comp.S~=0
                fun_block= [fun_block (arrayfun(@(i) trace((data.F_s(:,:,i)/X)*(data.G_s(:,:,i)/X)), 1:no_comp.S))]; 
            end
            if no_comp.M~=0
                fun_block= [fun_block (arrayfun(@(i) trace(data.P_m(:,:,i)*X*(data.Q_m(:,:,i)/X)), 1:no_comp.M))]; 
            end
            fun_block=fun_block';
%%
       f=user_fun(fun_block);
end