function[no_comp]=number_of_comp(data)
%% no of components
        if (isfield(data,'C_p'))
            no_comp.P=size(data.C_p,3);
        else
            no_comp.P=0;
        end
        if (isfield(data,'D_q'))
            no_comp.Q=size(data.D_q,3);
        else
            no_comp.Q=0;
        end
        %%%$$$$$%%%
        if data.logdt==1
            no_comp.logdt=1;
        elseif data.logdt==0
            no_comp.logdt=0;
        end
        %%%$$$$$%%%    
        if (isfield(data,'A_r'))
            no_comp.R=size(data.A_r,3);
        else
            no_comp.R=0;
        end
        if (isfield(data,'F_s'))
            no_comp.S=size(data.F_s,3);
        else
            no_comp.S=0;
        end
        if (isfield(data,'P_m'))
            no_comp.M=size(data.P_m,3);
        else
            no_comp.M=0;

        end
        %