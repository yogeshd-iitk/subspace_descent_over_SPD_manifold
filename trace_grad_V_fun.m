function[trace_grad_V]=trace_grad_V_fun(h,Int_var,dir_idx)
%% for offdiagonal directions
    if (isfield(dir_idx,'offd_col')) 
        trace_grad_V.offd=Int_var.hess_dh_offd*h;
    end
        
 %% for diagonal direction 
    if (isfield(dir_idx,'diag_row'))
        trace_grad_V.diag=Int_var.hess_dh_diag*h;
    end
end