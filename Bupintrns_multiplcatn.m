  function [M]=Bupintrns_multiplcatn(W,b_offd,d,dir_idx,M)      
        F_1 =zeros(W,W,size(M,3));
        if (isfield(dir_idx,'offd_col'))
            F_1(:,dir_idx.offd_row,:)=-M(:,dir_idx.offd_col,:).*(b_offd');
            M(:,[dir_idx.offd_col dir_idx.offd_row],:)=M(:,[dir_idx.offd_col dir_idx.offd_row],:)./d(1,[dir_idx.offd_col dir_idx.offd_row]);
        end
        if (isfield(dir_idx,'diag_row'))
            M(:,[dir_idx.diag_row],:)=M(:,[dir_idx.diag_row],:)./d(1,[dir_idx.diag_row]);
        else
        end
            M=M+F_1; 
  end