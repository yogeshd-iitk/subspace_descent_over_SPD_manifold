%%  Greedy subspace selection as specified in the paper 
%
function [dir_idx]=subspace_select_second_order_greedy(W,Intrm_var_mat,no_comp,h) 
a=[1:1:W^2];
A=reshape(a,[W,W]);
ind_offd=A(logical(tril(ones(W),-1)));
sz=[W,W];
[dir_idx_geedy.offd_row,dir_idx_geedy.offd_col] = ind2sub(sz,ind_offd);
dir_idx_geedy.diag_row=[1:W];
ind_diag=diag(A);
%%  Intermediate_variables calculation
    [Int_var]=intrmdt_var(dir_idx_geedy,Intrm_var_mat,no_comp);
%
%% beta calculation
    [trace_grad_V]=trace_grad_V_fun(h,Int_var,dir_idx_geedy);
    comp_para=zeros(W,W);  % comp_para is comparison parameter.: beta_ij 
    %
        if (isfield(dir_idx_geedy,'offd_col')) 
            comp_para_SD_offd=trace_grad_V.offd;
        end
        if (isfield(dir_idx_geedy,'diag_row'))
            comp_para_SD_diag=trace_grad_V.diag;
        end
    %
    comp_para(ind_offd)=comp_para_SD_offd;
    comp_para(ind_diag)=comp_para_SD_diag;
%
%% greedy subspace directions selection
AA=comp_para(:); % convert lower triangular matrix into column vector
nz_loc = tril(ones(W));  % location of entries of lower triangular part
nz_loc_vec=nz_loc(:);
AA_nz=nz_loc_vec>0; % logical vector
AA_ext=[AA, [1:1:W^2]'];
[~,idx] = sort(AA_ext(:,1),'ComparisonMethod','abs');
AA_ext_sorted=AA_ext(idx,:);
AA_nz_sorted=AA_nz(idx);
[~,idx_rvrs] = sort(AA_ext_sorted(:,2),'ComparisonMethod','abs');
sz=[W W];
i=1;
k=0;
while k<W
        temp_ind=find(AA_nz_sorted,1,'last'); % selceting first maximum entry which does not lie on  already chosen row or column of the chosen entries
        ind_current=AA_ext_sorted(temp_ind,2);
%         i=i+1;
        [row,col] = ind2sub(sz,ind_current); % location of maximum element
        if row ==col
            IND(i)=ind_current;
            i=i+1;
            k=k+1;
        else
            idx_2=(comp_para(row,row)^2+comp_para(col,col)^2 < sqrt(2)*comp_para(row,col)^2); % value of off-diagonal is larger than the diagonal entries combined
            if idx_2
                IND(i)=ind_current;
                i=i+1;
            else
                IND(i:i+1)=[(col-1)*W+col (row-1)*W+row];
                i=i+2;
            end
            k=k+2;
        end
        Ind=[col:W:col*W row:W:row*W, (col-1)*W+col:1:col*W, (row-1)*W+row+1:1:row*W]; % collect the other indices falling into the same row and column of the chosen entry
        Ind(col)=[]; % remove repeated  entry corresponding to chosen entry
        AA_nz_sorted(idx_rvrs(Ind))=false;
end
%
IND=sort(IND);
[final(:,2),final(:,1)]= ind2sub(sz,IND);
diag_idx=final(:,1)==final(:,2);
r_7=final(diag_idx,1)';
R_3=final(~diag_idx,:)'; % R_3 is already sorted according to first row as IND is sorted
%% storing indices
dir_idx.offd_col=R_3(1,:);
dir_idx.offd_row=R_3(2,:);
dir_idx.diag_row=r_7;     
end
