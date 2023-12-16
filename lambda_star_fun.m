%% lambda_star values corresponding to chosen subspace directions
function [lambda_star]=lambda_star_fun(dir_idx,Intrm_var_mat,no_comp,h,fun_block)
%%  Intermediate_variables calculation
        [Int_var]=intrmdt_var(dir_idx,Intrm_var_mat,no_comp);
    %
    %% directional derivative calculation
        [dir_deri]=direct_deriv(Int_var,dir_idx,fun_block);
    %
    %% alpha_SD calculation
        [trace_grad_V]=trace_grad_V_fun(h,Int_var,dir_idx);
        %
        [trace_V_grad_V_Xin]=trace_V_grad_V_Xin_fun(h,Int_var,dir_idx);
        %
        [trace_hessV_V]=trace_hessV_V_fun(h,dir_deri,Int_var,dir_idx,no_comp);
        %
       if (isfield(dir_idx,'offd_col')) 
            lambda_star.offd=-(trace_grad_V.offd)./(trace_hessV_V.offd+trace_V_grad_V_Xin.offd);
       end
    %
       if (isfield(dir_idx,'diag_row'))
            lambda_star.diag=-(trace_grad_V.diag)./(trace_hessV_V.diag+trace_V_grad_V_Xin.diag);
       end