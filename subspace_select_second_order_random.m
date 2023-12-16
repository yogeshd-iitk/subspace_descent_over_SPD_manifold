%% Subspace selection for RRSD (as specified in the paper in the section "Convergence analysis")
function [dir_idx]=subspace_select_second_order_random(W,subspace_dimension) 
    if mod(W,2)==0
        t=0;
    else
        t=(W+1)/2;
    end
%
    if strcmp(subspace_dimension, 'multi')  % multidimensional subspace
            p=1/W;
            if mod(W,2)==0
                    r_1=randperm(W);    % random permutation of the set {1,2,...W}  for random selection of Eij
                    r_3=reshape(r_1,2,[]); % selecting the directions
                    r_3=sort(r_3);
                    r_2=rand(W/2,1); 
                    r_2=(r_2<p); % indicator to select between diagonal or off diagonal direction
                    idx_3= r_2>0;
                    r_5=r_3(:,idx_3);
                    dir_idx.diag_row=[r_5(1,:) r_5(2,:)];
                    R_2=~r_2;
                    R_3=r_3(:,R_2>0);
                    [~,idx1_4]=sort(R_3(1,:));
                    R_3=R_3(:,idx1_4);
                    dir_idx.offd_col=R_3(1,:);
                    dir_idx.offd_row=R_3(2,:);
            else
                    r_1=randperm(W);
                    r_5=r_1(:,[1:W-1]);
                    r_3=reshape(r_5,2,[]);
                    r_3(:,(W+1)/2)=r_1(1,W)*ones(2,1);
                    r_3=sort(r_3);
                    r_2=rand((W-1)/2,1);
                    r_2=(r_2<p);
                    r_2((W+1)/2)=1;
                    r_6= r_3(:,1:end-1);
                    r_8=r_2(1:end-1,:);
                    idx_3= r_8>0;
                    r_5=r_6(:,idx_3);
                    dir_idx.diag_row=[r_5(1,:) r_5(2,:) r_3(1,t)];
                    R_2=~r_2;
                    R_3=r_3(:,R_2>0);
                    [~,idx1_4]=sort(R_3(1,:));
                    R_3=R_3(:,idx1_4);
                    dir_idx.offd_col=R_3(1,:);
                    dir_idx.offd_row=R_3(2,:);
            end
    %
    elseif  strcmp(subspace_dimension, 'uni')   % one-dimensional subspace
                R_1=randi(W);
                R_2=randi(W);
                if R_1~=R_2
                    dir_idx.offd_col=min(R_1,R_2); %column
                    dir_idx.offd_row=max(R_1,R_2); % row
                else
                    dir_idx.diag_row=R_1;
                end
    end
end