function[Intrm_var_mat,f,fun_block,fun_temp]=inti_var(B,Bin,latest,identity_initi,no_comp,data)
fun_temp=[];
%% For identity initialization
        if identity_initi==1
            if no_comp.P~=0
                Intrm_var_mat.M_1p=data.C_p;
            end
            if no_comp.Q~=0
                Intrm_var_mat.M_2q=data.D_q;
            end
            if no_comp.R~=0
                Intrm_var_mat.M_41r=data.A_r;
                Intrm_var_mat.M_42r=data.H_r;
            end
            if no_comp.S~=0
                Intrm_var_mat.M_51s=data.F_s;
                Intrm_var_mat.M_52s=data.G_s;
            end
            if no_comp.M~=0
                Intrm_var_mat.M_61m=data.P_m;
                Intrm_var_mat.M_62m=data.Q_m;
            end
        else
            %% for MATLAB R2020b or above
            if latest==1
                if no_comp.P~=0
                    Intrm_var_mat.M_1p=Pagemtimes(Pagemtimes(Bin,data.C_p),Bin');
                end
                if no_comp.Q~=0
                    Intrm_var_mat.M_2q=Pagemtimes(Pagemtimes(B',data.D_q),B);
                end
                if no_comp.R~=0
                    Intrm_var_mat.M_41r=Pagemtimes(Pagemtimes(B',data.A_r),B);
                    Intrm_var_mat.M_42r=Pagemtimes(Pagemtimes(B',data.H_r),B);
                end

                if no_comp.S~=0
                    Intrm_var_mat.M_51s=Pagemtimes(Pagemtimes(Bin,data.F_s),Bin');
                    Intrm_var_mat.M_52s=Pagemtimes(Pagemtimes(Bin,data.G_s),Bin');
                end

                if no_comp.M~=0
                    Intrm_var_mat.M_61m=Pagemtimes(Pagemtimes(Bin,data.P_m),B);
                    Intrm_var_mat.M_62m=Pagemtimes(Pagemtimes(Bin,data.Q_m),B);
                end
            else
            %%  for earlier versions of MATLAB where pagetimes doesnot work
                %
                for p=1:no_comp.P 
                    Intrm_var_mat.M_1p(:,:,p)= Bin*data.C_p(:,:,p)*(Bin'); %M_2= Bin*C*(Bin'); 
                end
                %
                for q=1:no_comp.Q 
                    Intrm_var_mat.M_2q(:,:,q)=(B')*data.D_q(:,:,q)*B; 
                end
                %
                for r=1:no_comp.R 
                    Intrm_var_mat.M_41r(:,:,r)=(B')*data.A_r(:,:,r)*B;
                    Intrm_var_mat.M_42r(:,:,r)=(B')*data.H_r(:,:,r)*B;
                end
                %
                for s=1:no_comp.S
                    Intrm_var_mat.M_51s(:,:,s)= Bin*data.F_s(:,:,s)*(Bin');
                    Intrm_var_mat.M_52s(:,:,s)= Bin*data.G_s(:,:,s)*(Bin');
                end
                %
                for m=1:no_comp.M 
                    Intrm_var_mat.M_61m(:,:,m)=(Bin*data.P_m(:,:,m))*B;
                    Intrm_var_mat.M_62m(:,:,m)=(Bin*data.Q_m(:,:,m))*B;
                end

            end
        end

        %% g function block calculations   
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
                fun_temp.uni_4r=arrayfun(@(i) trace(Intrm_var_mat.M_41r(:,:,i)*Intrm_var_mat.M_42r(:,:,i)), 1:no_comp.R);
                fun_block= [fun_block fun_temp.uni_4r]; 
            end
            if no_comp.S~=0
                fun_temp.uni_5s=arrayfun(@(i) trace(Intrm_var_mat.M_51s(:,:,i)*Intrm_var_mat.M_52s(:,:,i)), 1:no_comp.S);
                fun_block= [fun_block fun_temp.uni_5s];
            end
            if no_comp.M~=0
                fun_temp.uni_6m=arrayfun(@(i) trace(Intrm_var_mat.M_61m(:,:,i)*(Intrm_var_mat.M_62m(:,:,i)')), 1:no_comp.M);
                fun_block= [fun_block fun_temp.uni_6m];
            end
            fun_block=fun_block';
            % input('press enter')
    %% Function evaluation
        f=user_fun(fun_block);
 end
