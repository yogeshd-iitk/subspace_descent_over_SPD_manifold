function[trace_hessV_V]=trace_hessV_V_fun(h,dir_deri,Int_var,dir_idx,no_comp)
    if no_comp.Q~=0
            h(no_comp.P+1:no_comp.P+no_comp.Q)=[];
   end
%% for offdiagonal directions
    if (isfield(dir_idx,'offd_col')) 
       trace_hessV_V.offd=sum(dir_deri.Dh_offd.*Int_var.hess_dh_offd,2)+Int_var.hess_h_offd*h;
    end
        
 %% for diagonal directions 
    if (isfield(dir_idx,'diag_row'))
       trace_hessV_V.diag=sum(dir_deri.Dh_diag.*Int_var.hess_dh_diag,2)+Int_var.hess_h_diag*h;
    end
end