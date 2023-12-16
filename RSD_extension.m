%% Code for multi RRSD, uni RRSD and multi RGSD. 
% note that uni RGSD is not implemented in this code.
clc;
clear variables;
close all;
format long
W=input('enter the dimension of matrices='); 
iter=input('enter the no of iteration='); 
stepsz_method= 'adaptive'; % input('stepsize method='); %  'constant' ; % 
if  strcmp(stepsz_method, 'constant')
    stepsize= input('enter the value of constant stepsize=');
end
scale_step=1;  % scale the stepsize
latest=0 ; % 1 if MATLAB  R2020b or above is available
paper=0;    % for the seetings used in paper, optimum value is known, set paper =1 to plot results relative to the optimum value 
subspace_dimension='multi';    %  'uni';  %          input('subspace_dimension='); % 
algo= 'greedy';  % 'randomized'; %       input('algo='); %  
%% Initialization
    X_0=eye(W);
    if X_0==eye(W)
        identity_initi=1;
        B=eye(W); % chol(X_0)';    %%%%%%%%%%%%%%% B= identity matrix (I) when X_0= I
        Bin=eye(W);                %%%%%%%%%%%%%%% Bin= identity matrix (I) when X_0= I
    else
        identity_initi=0;
        % define B matrix
        % define Bin matrix
    end
%% data given by user
    %
    [data,no_comp,f_user]=user_input(W);
    %
%% Save data
    filename = ['user_data.mat'];
    save(filename,'W','X_0','B','Bin','data','no_comp','identity_initi');

%% symbolic gradient and Hessian
gradf = jacobian(f_user,symvar(f_user)).';
hessf = jacobian(gradf,symvar(f_user));
%
%% symbolic function form  to function file conversion
%
matlabFunction(f_user,"File","user_fun","Vars",{symvar(f_user).'},"Optimize",false);
%
%% symbolic gradient vector form  to function file conversion
%
matlabFunction(gradf,"File","user_fun_grad","Vars",{symvar(f_user).'},"Optimize",false);
%
%% symbolic Hessian matrix form  to function file conversion
%
matlabFunction(hessf,"File","user_fun_hess","Vars",{symvar(f_user).'},"Optimize",false);
%
%% Initialization 
result.f_second_order=zeros(iter+1,1);
result.T_second_order=zeros(iter,1);
%
%% intermediate variables initialization
[Intrm_var_mat,result.f_second_order(1),fun_block,fun_temp]=inti_var(B,Bin,latest,identity_initi,no_comp,data);
identity_initi=0;
Bin=[];
%
%% main algo
for m_2=1:iter
        
    tic
        m_2

    %% gradient evaluation w.r.t function blocks as variables
        h=user_fun_grad(fun_block);
    %
    %% subspace selection
            if strcmp(algo, 'randomized')
                [dir_idx]=subspace_select_second_order_random(W,subspace_dimension);
                if (isfield(dir_idx,'offd_col'))
                    result.offd_dir_random(m_2)=length(dir_idx.offd_col);   % counting no. of offdiagonal direction
                else
                    result.offd_dir_random(m_2)=0;
                end
                if (isfield(dir_idx,'diag_row'))
                    result.diag_dir_random(m_2)=length(dir_idx.diag_row);   % counting no. of  diagonal direction
                else
                    result.diag_dir_random(m_2)=0;
                end
            elseif strcmp(algo, 'greedy')
                [dir_idx]=subspace_select_second_order_greedy(W,Intrm_var_mat,no_comp,h);
                if (isfield(dir_idx,'offd_col'))
                    result.offd_dir_greedy(m_2)=length(dir_idx.offd_col);   % counting no. of offdiagonal direction
                else
                    result.offd_dir_greedy(m_2)=0;
                end
                if (isfield(dir_idx,'diag_row'))
                    result.diag_dir_greedy(m_2)=length(dir_idx.diag_row);   % counting no. of  diagonal direction
                else
                    result.diag_dir_greedy(m_2)=0;
                end
            end
    % step-size 
    if strcmp(stepsz_method, 'adaptive')
        [lambda_star]=lambda_star_fun(dir_idx,Intrm_var_mat,no_comp,h,fun_block);
            if (isfield(dir_idx,'offd_col')) % off_diagoal directions
                alpha_SD_offd=(1/scale_step)*lambda_star.offd;
            end
            if (isfield(dir_idx,'diag_row')) % Diagonal directions
                alpha_SD_diag=(1/scale_step)*lambda_star.diag;
            end
            %
    elseif strcmp(stepsz_method, 'constant')
            if (isfield(dir_idx,'offd_col')) % off_diagoal directions
                alpha_SD_offd=stepsize*ones(size(dir_idx.offd_col),1);
            end
            if (isfield(dir_idx,'diag_row')) % Diagonal directions
                alpha_SD_diag=stepsize*ones(size(dir_idx.diag_row),1);
            end
    end
    %% calculating update matrix B_{t+1}^{up}
        d=ones(1,W);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (isfield(dir_idx,'offd_col'))
            a=exp((1/sqrt(2))*alpha_SD_offd); 
            alpha_SD_offd=[];
            b=1./a;
            t_1=sqrt((a+b)/2);
            b_offd=(a-b)./(2*t_1); % offdiagonal entries of Bup matrix 
            a=[];
            b=[];
            d([dir_idx.offd_col dir_idx.offd_row])=[t_1 1./t_1];
            t_1=[];
        else
            b_offd=0;
        end
        if (isfield(dir_idx,'diag_row'))
            d(dir_idx.diag_row)=exp(0.5*alpha_SD_diag); %diagonal entries of Bup matrix
            alpha_SD_diag=[];
        end
        %
%% updating intermediate matrices
%
    %% updating variables corresponding to C_p
        if no_comp.P~=0
        %% update Intrm_var_mat.M_1p  : [Bup_inv*(M_1p*Bup_inv')]
            [temp_1]=Bupintrns_multiplcatn(W,b_offd,d,dir_idx,Intrm_var_mat.M_1p);
            [Intrm_var_mat.M_1p]=permute((Bupintrns_multiplcatn(W,b_offd,d,dir_idx,permute(temp_1,[2 1 3]))),[2 1 3]); 
            temp_1=[];
            %
        end
    %% updating variables corresponding to D_q
        if no_comp.Q~=0
        %% update Intrm_var_mat.M_2q : [Bup'*(M_1*Bup)=(((M_1*Bup)')*Bup)']
            [temp_2]=Bup_multiplcatn(W,b_offd,d,dir_idx,Intrm_var_mat.M_2q);
            [Intrm_var_mat.M_2q]=permute((Bup_multiplcatn(W,b_offd,d,dir_idx,permute(temp_2,[2 1 3]))),[2 1 3]);
            temp_2=[];
            %
        end
    %% updating variables corresponding to A_r and H_r
        if no_comp.R~=0
        if strcmp(subspace_dimension,'uni')
            [mod_entry_old.uni_41r]=mod_entry(dir_idx,Intrm_var_mat.M_41r);
            [mod_entry_old.uni_42r]=mod_entry(dir_idx,Intrm_var_mat.M_42r);
        end
        %% update Intrm_var_mat.M_41r : [Bup'*(M_41r*Bup)=(((M_1*Bup)')*Bup)']
            [temp_41]=Bup_multiplcatn(W,b_offd,d,dir_idx,Intrm_var_mat.M_41r);
            [Intrm_var_mat.M_41r]=permute((Bup_multiplcatn(W,b_offd,d,dir_idx,permute(temp_41,[2 1 3]))),[2 1 3]);
            temp_41=[];
            %
        %% update Intrm_var_mat.M_42r : [Bup'*(M_42r*Bup)=(((M_1*Bup)')*Bup)']
            [temp_42]=Bup_multiplcatn(W,b_offd,d,dir_idx,Intrm_var_mat.M_42r);
            [Intrm_var_mat.M_42r]=permute((Bup_multiplcatn(W,b_offd,d,dir_idx,permute(temp_42,[2 1 3]))),[2 1 3]);
            temp_42=[];
            %
        if strcmp(subspace_dimension,'uni')
            [mod_entry_new.uni_41r]=mod_entry(dir_idx,Intrm_var_mat.M_41r);
            [mod_entry_new.uni_42r]=mod_entry(dir_idx,Intrm_var_mat.M_42r);
        end
        end
    %% updating variables corresponding to F_s and G_s
        if no_comp.S~=0
        if strcmp(subspace_dimension,'uni')
            [mod_entry_old.uni_51s]=mod_entry(dir_idx,Intrm_var_mat.M_51s);
            [mod_entry_old.uni_52s]=mod_entry(dir_idx,Intrm_var_mat.M_52s);
        end
        %% update Intrm_var_mat.M_51s : [Bup_inv*(M_51s*Bup_inv')]
            [temp_51]=Bupintrns_multiplcatn(W,b_offd,d,dir_idx,Intrm_var_mat.M_51s);
            [Intrm_var_mat.M_51s]=permute((Bupintrns_multiplcatn(W,b_offd,d,dir_idx,permute(temp_51,[2 1 3]))),[2 1 3]); 
            temp_51=[];
            %
        %% update Intrm_var_mat.M_52s : [Bup_inv*(M_52s*Bup_inv')]
            [temp_52]=Bupintrns_multiplcatn(W,b_offd,d,dir_idx,Intrm_var_mat.M_52s);
            [Intrm_var_mat.M_52s]=permute((Bupintrns_multiplcatn(W,b_offd,d,dir_idx,permute(temp_52,[2 1 3]))),[2 1 3]); 
            temp_52=[];
            %
        if strcmp(subspace_dimension,'uni')
            [mod_entry_new.uni_51s]=mod_entry(dir_idx,Intrm_var_mat.M_51s);
            [mod_entry_new.uni_52s]=mod_entry(dir_idx,Intrm_var_mat.M_52s);
        end
        end
    %% updating variables corresponding to P_m and Q_m
        if no_comp.M~=0
        if strcmp(subspace_dimension,'uni')
            [mod_entry_old.uni_61m]=mod_entry(dir_idx,Intrm_var_mat.M_61m);
            [mod_entry_old.uni_62m]=mod_entry(dir_idx,Intrm_var_mat.M_62m);
        end
            %% update Intrm_var_mat.M_61m : [Bup_inv*(M_61m*Bup)]
            [temp_61]=Bup_multiplcatn(W,b_offd,d,dir_idx,Intrm_var_mat.M_61m);
            [Intrm_var_mat.M_61m]=permute((Bupintrns_multiplcatn(W,b_offd,d,dir_idx,permute(temp_61,[2 1 3]))),[2 1 3]); 
            temp_61=[];
            %
            %% update Intrm_var_mat.M_62m : [Bup_inv*(M_62m*Bup)]  
            [temp_62]=Bup_multiplcatn(W,b_offd,d,dir_idx,Intrm_var_mat.M_62m);
            [Intrm_var_mat.M_62m]=permute((Bupintrns_multiplcatn(W,b_offd,d,dir_idx,permute(temp_62,[2 1 3]))),[2 1 3]); 
            temp_62=[];
            %
        if strcmp(subspace_dimension,'uni')
            [mod_entry_new.uni_61m]=mod_entry(dir_idx,Intrm_var_mat.M_61m);
            [mod_entry_new.uni_62m]=mod_entry(dir_idx,Intrm_var_mat.M_62m);
        end
        end
        %% update B matrix (B=B*Bup)
        B=Bup_multiplcatn(W,b_offd,d,dir_idx,B);
        d=[];
    %%  g function block calculations   
    % numell=sum(cell2mat(struct2cell(no_comp),'all');
    % fun_block=zeros(1,numell);
            if no_comp.P~=0
                fun_block=arrayfun(@(i) trace(Intrm_var_mat.M_1p(:,:,i)), 1:no_comp.P); % g_1p % g_1P=trace(Intrm_var_mat.M_1p);
            else
                fun_block=[];
            end
            if no_comp.Q~=0
                fun_block=[fun_block (arrayfun(@(i)  trace(Intrm_var_mat.M_2q(:,:,i)), 1:no_comp.Q))];%trace(Intrm_var_mat.M_2q);
            end
            if no_comp.logdt~=0
                fun_block= [fun_block 2*sum(log(diag(B)))]; 
            end
                if no_comp.R~=0
                    if strcmp(subspace_dimension,'uni')
                        fun_temp.uni_4r= fun_temp.uni_4r +sum((mod_entry_new.uni_41r.*mod_entry_new.uni_42r),1)-sum((mod_entry_old.uni_41r.*mod_entry_old.uni_42r),1);
                    else
                        fun_temp.uni_4r=arrayfun(@(i) trace(Intrm_var_mat.M_41r(:,:,i)*Intrm_var_mat.M_42r(:,:,i)), 1:no_comp.R);
                    end
                    fun_block= [fun_block (fun_temp.uni_4r)]; % trace(Intrm_var_mat.M_41r*Intrm_var_mat.M_42r);
                end
                if no_comp.S~=0
                    if strcmp(subspace_dimension,'uni')
                        fun_temp.uni_5s= fun_temp.uni_5s +sum((mod_entry_new.uni_51s.*mod_entry_new.uni_52s),1)-sum((mod_entry_old.uni_51s.*mod_entry_old.uni_52s),1);
                    else
                        fun_temp.uni_5s=arrayfun(@(i) trace(Intrm_var_mat.M_51s(:,:,i)*Intrm_var_mat.M_52s(:,:,i)), 1:no_comp.S);
                    end
                    fun_block= [fun_block (fun_temp.uni_5s)]; % trace(Intrm_var_mat.M_51s*Intrm_var_mat.M_52s);
                end
                if no_comp.M~=0
                    if strcmp(subspace_dimension,'uni')
                        fun_temp.uni_6m= fun_temp.uni_6m +sum((mod_entry_new.uni_61m.*mod_entry_new.uni_62m),1)-sum((mod_entry_old.uni_61m.*mod_entry_old.uni_62m),1);
                    else
                        fun_temp.uni_6m=arrayfun(@(i) trace(Intrm_var_mat.M_61m(:,:,i)*(Intrm_var_mat.M_62m(:,:,i)')), 1:no_comp.M);
                    end
                    fun_block= [fun_block (fun_temp.uni_6m)]; % trace(Intrm_var_mat.M_61m*Intrm_var_mat.M_62m);
                end
            fun_block=fun_block';
            % input('press enter')
    %% Function evaluation
        result.f_second_order(m_2+1)=user_fun(fun_block);
        %
        f_iter=result.f_second_order(m_2+1)
        toc
        result.T_second_order(m_2,1)=toc; 
        clc;
end
%% Saving and plotting results for setting used in paper
if paper==1
% Evaluate optimal function value for the setting  used in paper
    k=sign(gradf(no_comp.P+no_comp.Q+1));
    if k==0
        X_star=sqrtm(data.C_p);
        f_x_star=trace(X_star\data.C_p)+trace(data.D_q*X_star)
    elseif k==-1
        [V,Egn] = eig(data.C_p);
        egn_c=diag(Egn);
        egn_x=(1+sqrt(1+4*egn_c))/2;
        X_star=V*diag(egn_x)*(V');
        f_x_star=trace(X_star\data.C_p)+trace(data.D_q*X_star)-log(det(X_star))
    end
% save results 
        if strcmp(algo, 'randomized')
            filename = ['function_value_second_order_' num2str(W) '_' num2str(iter) '_' num2str(scale_step)  '_', subspace_dimension,'RRSD' '.mat'];
        elseif strcmp(algo, 'greedy') 
            filename = ['function_value_second_order_' num2str(W) '_' num2str(iter) '_' num2str(scale_step) '_', subspace_dimension,'RGSD' '.mat'];
        end
            save(filename,'W','X_0','B','data','result','-v7.3');    
%  plot results
            load(filename,'result'); 
             f_result=result.f_second_order;
             plot([0:iter]',f_result-f_x_star,'b-*')
             ylabel({'$f(X)-f(X^{\star})$'},'Interpreter','latex')
             xlabel('Iteration')
end
%% Plot results 
if paper==0
% plotting function values
     figure
     f_result=result.f_second_order;
     plot([0:iter]',f_result,'b-*')
     ylabel({'$f(X)$'},'Interpreter','latex')
     xlabel('Iteration')
 %% save results 
    if strcmp(algo, 'randomized')
        filename = ['function_value_second_order_' num2str(W) '_' num2str(iter) '_' num2str(scale_step)  '_', subspace_dimension,'RRSD' '.mat'];
    elseif strcmp(algo, 'greedy') 
        filename = ['function_value_second_order_' num2str(W) '_' num2str(iter) '_' num2str(scale_step) '_', subspace_dimension,'RGSD' '.mat'];
    end
        save(filename,'W','X_0','B','data','result','-v7.3');   
 %%
end
 result.f_second_order(end)
sum(result.T_second_order)
