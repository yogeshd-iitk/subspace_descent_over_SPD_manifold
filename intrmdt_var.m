function [Int_var]=intrmdt_var(dir_idx,Intrm_var_mat,no_comp)
%% offdiagonal terms
    if (isfield(dir_idx,'offd_col'))
        l_offd=length(dir_idx.offd_col); 
        if no_comp.P~=0
            for p=1:no_comp.P
                T_M_1p=Intrm_var_mat.M_1p(:,:,p); %%%%%%%%%%%%%%%%%%%%%%%%
                T_1_offd(:,p,1)=(sqrt(2)*T_M_1p(sub2ind(size(T_M_1p),dir_idx.offd_row,dir_idx.offd_col)))'; % T_1_offd(:,p,1)=T_11_offd(:,p)
                T_1_offd(:,p,2)=((1/2)*(T_M_1p(sub2ind(size(T_M_1p),dir_idx.offd_row,dir_idx.offd_row))+T_M_1p(sub2ind(size(T_M_1p),dir_idx.offd_col,dir_idx.offd_col))))'; % T_1_offd(:,p,2)=T_12_offd(:,p)
            end
            Int_var.hess_dh_offd=-T_1_offd(:,:,1);
            Int_var.hess_h_offd=2*T_1_offd(:,:,2);
            Int_var.grad_offd=-T_1_offd(:,:,2);
        else
            Int_var.hess_dh_offd=[];
            Int_var.hess_h_offd=[];
            Int_var.grad_offd=[];
        end
        if no_comp.Q~=0
            for q=1:no_comp.Q
                T_M_2q=Intrm_var_mat.M_2q(:,:,q);
                T_2_offd(:,q,1)=(sqrt(2)* T_M_2q(sub2ind(size(T_M_2q),dir_idx.offd_row,dir_idx.offd_col)))'; %  T_2_offd(:,q,1)=T_21_offd(:,q)
                T_2_offd(:,q,2)=((1/2)*(T_M_2q(sub2ind(size(T_M_2q),dir_idx.offd_row,dir_idx.offd_row))+T_M_2q(sub2ind(size(T_M_2q),dir_idx.offd_col,dir_idx.offd_col))))'; % T_2_offd(:,q,2)=T_22_offd(:,q)
            end
            Int_var.hess_dh_offd=[Int_var.hess_dh_offd  T_2_offd(:,:,1)];
            Int_var.grad_offd=[Int_var.grad_offd  T_2_offd(:,:,2)];
        end
            if no_comp.logdt~=0
                T_31_offd=zeros(l_offd,1); %%%%%%%%%%%%%%%%%%%%%%%%
                T_32_offd=ones(l_offd,1);%%%%%%%%%%%%%%%%%%%%%%%%
                Int_var.hess_dh_offd=[Int_var.hess_dh_offd  T_31_offd];
                Int_var.hess_h_offd=[Int_var.hess_h_offd  (-T_32_offd)];
                Int_var.grad_offd=[Int_var.grad_offd  T_32_offd];
            end
        if no_comp.R~=0
            for r=1:no_comp.R
                T_M_41r=Intrm_var_mat.M_41r(:,:,r);
                T_M_42r=Intrm_var_mat.M_42r(:,:,r);
                T_4_offd(:,r,1)=(1/sqrt(2))*sum((T_M_41r(dir_idx.offd_col,:).*T_M_42r(dir_idx.offd_row,:)+T_M_42r(dir_idx.offd_col,:).*T_M_41r(dir_idx.offd_row,:)),2); % T_4_offd(:,r,1)=T_41_offd(:,r)
                T_4_offd(:,r,2)=((1/2)*(2*T_M_41r(sub2ind(size(T_M_41r),dir_idx.offd_row,dir_idx.offd_col)).*T_M_42r(sub2ind(size(T_M_42r),dir_idx.offd_row,dir_idx.offd_col))+T_M_41r(sub2ind(size(T_M_41r),dir_idx.offd_row,dir_idx.offd_row)).*T_M_42r(sub2ind(size(T_M_42r),dir_idx.offd_col,dir_idx.offd_col))+T_M_41r(sub2ind(size(T_M_41r),dir_idx.offd_col,dir_idx.offd_col)).*T_M_42r(sub2ind(size(T_M_42r),dir_idx.offd_row,dir_idx.offd_row))))';% T_4_offd(:,r,2)=T_42_offd(:,r)
                T_4_offd(:,r,3)=(1/2)*sum((T_M_41r(dir_idx.offd_row,:).*T_M_42r(dir_idx.offd_row,:)+T_M_42r(dir_idx.offd_col,:).*T_M_41r(dir_idx.offd_col,:)),2); % T_4_offd(:,r,3)=T_43_offd(:,r)
            end
            Int_var.hess_dh_offd=[Int_var.hess_dh_offd  (2*T_4_offd(:,:,1))];
            Int_var.hess_h_offd=[Int_var.hess_h_offd  (2*T_4_offd(:,:,2))];
            Int_var.grad_offd=[Int_var.grad_offd  (2*T_4_offd(:,:,3))];
        end
        if no_comp.S~=0
            for s=1:no_comp.S
                T_M_51s=Intrm_var_mat.M_51s(:,:,s);
                T_M_52s=Intrm_var_mat.M_52s(:,:,s);
                T_5_offd(:,s,1)=(1/sqrt(2))*sum((T_M_51s(dir_idx.offd_col,:).*T_M_52s(dir_idx.offd_row,:)+T_M_52s(dir_idx.offd_col,:).*T_M_51s(dir_idx.offd_row,:)),2); % T_5_offd(:,s,1)=T_51_offd(:,s)
                T_5_offd(:,s,2)=(1/2)*sum((T_M_51s(dir_idx.offd_row,:).*T_M_52s(dir_idx.offd_row,:)+T_M_52s(dir_idx.offd_col,:).*T_M_51s(dir_idx.offd_col,:)),2)  ;  %T_5_offd(:,s,2)=T_52_offd(:,s)
                T_5_offd(:,s,3)=((1/2)*(2*T_M_51s(sub2ind(size(T_M_51s),dir_idx.offd_row,dir_idx.offd_col)).*T_M_52s(sub2ind(size(T_M_52s),dir_idx.offd_row,dir_idx.offd_col))+T_M_51s(sub2ind(size(T_M_51s),dir_idx.offd_row,dir_idx.offd_row)).*T_M_52s(sub2ind(size(T_M_52s),dir_idx.offd_col,dir_idx.offd_col))+T_M_51s(sub2ind(size(T_M_51s),dir_idx.offd_col,dir_idx.offd_col)).*T_M_52s(sub2ind(size(T_M_52s),dir_idx.offd_row,dir_idx.offd_row))))'; %T_5_offd(:,s,3)=T_53_offd(:,s)
            end
            Int_var.hess_dh_offd=[Int_var.hess_dh_offd (-2*T_5_offd(:,:,1))];
            Int_var.hess_h_offd=[Int_var.hess_h_offd  (4*T_5_offd(:,:,2)+2*T_5_offd(:,:,3))];
            Int_var.grad_offd=[Int_var.grad_offd  (-2*T_5_offd(:,:,2))];
        end
        if no_comp.M~=0
            for m=1:no_comp.M
                T_M_61m=Intrm_var_mat.M_61m(:,:,m);
                T_M_62m=Intrm_var_mat.M_62m(:,:,m);
                T_6_offd(:,m,1)=(1/sqrt(2))*sum((T_M_62m(:,dir_idx.offd_col).*T_M_61m(:,dir_idx.offd_row)+T_M_62m(:,dir_idx.offd_row).*T_M_61m(:,dir_idx.offd_col)),1)'; %T_6_offd(:,m,1)=T_61_offd(:,m)
                T_6_offd(:,m,2)=(1/sqrt(2))*sum((T_M_61m(dir_idx.offd_col,:).*T_M_62m(dir_idx.offd_row,:)+T_M_61m(dir_idx.offd_row,:).*T_M_62m(dir_idx.offd_col,:)),2); %T_6_offd(:,m,2)=T_62_offd(:,m)
                T_6_offd(:,m,3)=((1/2)*(T_M_61m(sub2ind(size(T_M_61m),dir_idx.offd_col,dir_idx.offd_row)).*T_M_62m(sub2ind(size(T_M_62m),dir_idx.offd_row,dir_idx.offd_col))+T_M_61m(sub2ind(size(T_M_61m),dir_idx.offd_row,dir_idx.offd_row)).*T_M_62m(sub2ind(size(T_M_62m),dir_idx.offd_col,dir_idx.offd_col))+T_M_61m(sub2ind(size(T_M_61m),dir_idx.offd_col,dir_idx.offd_col)).*T_M_62m(sub2ind(size(T_M_62m),dir_idx.offd_row,dir_idx.offd_row))+T_M_61m(sub2ind(size(T_M_61m),dir_idx.offd_row,dir_idx.offd_col)).*T_M_62m(sub2ind(size(T_M_62m),dir_idx.offd_col,dir_idx.offd_row))))'; % T_6_offd(:,m,3)=T_63_offd(:,m)
                T_6_offd(:,m,4)=(1/2)*sum((T_M_61m(dir_idx.offd_row,:).*T_M_62m(dir_idx.offd_row,:)+T_M_61m(dir_idx.offd_col,:).*T_M_62m(dir_idx.offd_col,:)),2); % T_6_offd(:,m,4)=T_64_offd(:,m)
                T_6_offd(:,m,5)=(1/2)*sum((T_M_62m(:,dir_idx.offd_row).*T_M_61m(:,dir_idx.offd_row)+T_M_62m(:,dir_idx.offd_col).*T_M_61m(:,dir_idx.offd_col)),1)'; % T_6_offd(:,m,5)=T_65_offd(:,m)
            end
            Int_var.hess_dh_offd=[Int_var.hess_dh_offd  (T_6_offd(:,:,1)-T_6_offd(:,:,2))];
            Int_var.hess_h_offd=[Int_var.hess_h_offd  (2*(T_6_offd(:,:,4)-T_6_offd(:,:,3)))];
            Int_var.grad_offd=[Int_var.grad_offd  (T_6_offd(:,:,5)-T_6_offd(:,:,4))];
        end
    end
%
%
%% diagonal terms 
    if (isfield(dir_idx,'diag_row'))
        l_diag=length(dir_idx.diag_row);  
        if no_comp.P~=0
            for p=1:no_comp.P
                T_M_1p=Intrm_var_mat.M_1p(:,:,p); %%%%%%%%%%%%%%%%%%%%%%%%
                T_1_diag(:,p,1)=(T_M_1p(sub2ind(size(T_M_1p),dir_idx.diag_row,dir_idx.diag_row)))';%T_1_offd(:,p,1)=T_11_diag(:,p)
                T_1_diag(:,p,2)=(T_M_1p(sub2ind(size(T_M_1p),dir_idx.diag_row,dir_idx.diag_row)))'; % T_1_offd(:,p,2)=T_12_diag(:,p)
            end
            Int_var.hess_dh_diag=-T_1_diag(:,:,1);
            Int_var.hess_h_diag=2*T_1_diag(:,:,2);
            Int_var.grad_diag=-T_1_diag(:,:,2);
        else
            Int_var.hess_dh_diag=[];
            Int_var.hess_h_diag=[];
            Int_var.grad_diag=[];
        end
        if no_comp.Q~=0
            for q=1:no_comp.Q
                T_M_2q=Intrm_var_mat.M_2q(:,:,q);
                T_2_diag(:,q,1)=(T_M_2q(sub2ind(size(T_M_2q),dir_idx.diag_row,dir_idx.diag_row)))'; % T_2_diag(:,q,1)=T_21_diag(:,q)
                T_2_diag(:,q,2)=(T_M_2q(sub2ind(size(T_M_2q),dir_idx.diag_row,dir_idx.diag_row)))'; %T_1_diag(:,q,2)=T_22_diag(:,q)
            end
            Int_var.hess_dh_diag=[Int_var.hess_dh_diag  T_2_diag(:,:,1)];
            Int_var.grad_diag=[Int_var.grad_diag  T_2_diag(:,:,2)];
        end
            if no_comp.logdt~=0
                T_31_diag=ones(l_diag,1); %%%$$$$$%%%
                T_32_diag=ones(l_diag,1); %%%$$$$$%%%
                Int_var.hess_dh_diag=[Int_var.hess_dh_diag  T_31_diag];
                Int_var.hess_h_diag=[Int_var.hess_h_diag  (-T_32_diag)];
                Int_var.grad_diag=[Int_var.grad_diag  T_32_diag];
            end
        if no_comp.R~=0
            for r=1:no_comp.R
                T_M_41r=Intrm_var_mat.M_41r(:,:,r);
                T_M_42r=Intrm_var_mat.M_42r(:,:,r);   
                T_4_diag(:,r,1)=sum((T_M_41r(dir_idx.diag_row,:).*T_M_42r(dir_idx.diag_row,:)),2); %T_4_diag(:,r,1)=T_41_diag(:,r)
                T_4_diag(:,r,2)=(T_M_41r(sub2ind(size(T_M_41r),dir_idx.diag_row,dir_idx.diag_row)).*T_M_42r(sub2ind(size(T_M_42r),dir_idx.diag_row,dir_idx.diag_row)))'; % T_4_diag(:,r,2)=T_42_diag(:,r)
                T_4_diag(:,r,3)=sum((T_M_41r(dir_idx.diag_row,:).*T_M_42r(dir_idx.diag_row,:)),2); %T_4_diag(:,r,3)=T_43_diag(:,r)
            end
            Int_var.hess_dh_diag=[Int_var.hess_dh_diag (2*T_4_diag(:,:,1))];
            Int_var.hess_h_diag=[Int_var.hess_h_diag  (2*T_4_diag(:,:,2))];
            Int_var.grad_diag=[Int_var.grad_diag  (2*T_4_diag(:,:,3))];
        end
        if no_comp.S~=0
            for s=1:no_comp.S
                T_M_51s=Intrm_var_mat.M_51s(:,:,s);
                T_M_52s=Intrm_var_mat.M_52s(:,:,s);   
                T_5_diag(:,s,1)=sum((T_M_51s(dir_idx.diag_row,:).*T_M_52s(dir_idx.diag_row,:)),2); % T_5_diag(:,s,1)=T_51_diag(:,s)
                T_5_diag(:,s,2)=sum((T_M_51s(dir_idx.diag_row,:).*T_M_52s(dir_idx.diag_row,:)),2); % T_5_diag(:,s,2)=T_52_diag(:,s)
                T_5_diag(:,s,3)=(T_M_51s(sub2ind(size(T_M_51s),dir_idx.diag_row,dir_idx.diag_row)).*T_M_52s(sub2ind(size(T_M_52s),dir_idx.diag_row,dir_idx.diag_row)))'; % T_5_diag(:,s,3)=T_53_diag(:,s)
            end
            Int_var.hess_dh_diag=[Int_var.hess_dh_diag (-2*T_5_diag(:,:,1))];
            Int_var.hess_h_diag=[Int_var.hess_h_diag  (4*T_5_diag(:,:,2)+2*T_5_diag(:,:,3))];
            Int_var.grad_diag=[Int_var.grad_diag  (-2*T_5_diag(:,:,2))];
        end
        if no_comp.M~=0
            for m=1:no_comp.M
                T_M_61m=Intrm_var_mat.M_61m(:,:,m);
                T_M_62m=Intrm_var_mat.M_62m(:,:,m);
                T_6_diag(:,m,1)=sum((T_M_62m(:,dir_idx.diag_row).*T_M_61m(:,dir_idx.diag_row)),1)'; % T_6_diag(:,m,1)=T_61_diag(:,m)
                T_6_diag(:,m,2)=sum((T_M_61m(dir_idx.diag_row,:).*T_M_62m(dir_idx.diag_row,:)),2); % T_6_diag(:,m,2)=T_62_diag(:,m)
                T_6_diag(:,m,3)=(T_M_61m(sub2ind(size(T_M_61m),dir_idx.diag_row,dir_idx.diag_row)).*T_M_62m(sub2ind(size(T_M_62m),dir_idx.diag_row,dir_idx.diag_row)))'; % T_6_diag(:,m,3)=T_63_diag(:,m)
                T_6_diag(:,m,4)=sum((T_M_61m(dir_idx.diag_row,:).*T_M_62m(dir_idx.diag_row,:)),2); % _6_diag(:,m,4)=T_64_diag(:,m)
                T_6_diag(:,m,5)=sum((T_M_62m(:,dir_idx.diag_row).*T_M_61m(:,dir_idx.diag_row)),1)'; %T_6_diag(:,m,5)=T_65_diag(:,m)
            end
            Int_var.hess_dh_diag=[Int_var.hess_dh_diag  (T_6_diag(:,:,1)-T_6_diag(:,:,2))];
            Int_var.hess_h_diag=[Int_var.hess_h_diag  (2*(T_6_diag(:,:,4)-T_6_diag(:,:,3)))];
            Int_var.grad_diag=[Int_var.grad_diag  (T_6_diag(:,:,5)-T_6_diag(:,:,4))];
        end
    end
%
%
end