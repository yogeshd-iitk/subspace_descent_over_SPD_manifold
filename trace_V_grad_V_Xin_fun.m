function[trace_V_grad_V_Xin]=trace_V_grad_V_Xin_fun(h,Int_var,dir_idx)
%% for offdiagonal directions
    if (isfield(dir_idx,'offd_col')) 
        trace_V_grad_V_Xin.offd=Int_var.grad_offd*h;
    end
        
 %% for diagonal direction   
     if (isfield(dir_idx,'diag_row'))
        trace_V_grad_V_Xin.diag=Int_var.grad_diag*h;
     end
end
