function[dir_deri]=direct_deriv(Int_var,dir_idx,fun_block)
%% directional derivative of h
hessf=user_fun_hess(fun_block);
% along off-diagonal direction
%
    if (isfield(dir_idx,'offd_col')) 
        dir_deri.Dh_offd= (hessf*(Int_var.hess_dh_offd'))';
    end
%
% along diagonal direction
%
    if (isfield(dir_idx,'diag_row'))
        dir_deri.Dh_diag= (hessf*(Int_var.hess_dh_diag'))';
    end
end