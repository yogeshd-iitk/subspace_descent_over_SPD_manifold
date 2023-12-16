function [mod_entry]=mod_entry(dir_idx,M)      
            if (isfield(dir_idx,'offd_col')) % off diagonal direction 
                mod_entry=squeeze(M(:,dir_idx.offd_col,:));
                mod_entry([dir_idx.offd_col,dir_idx.offd_row],:)=[];
                temp=squeeze(M(:,dir_idx.offd_row,:));
                temp([dir_idx.offd_col,dir_idx.offd_row],:)=[];
                if size(M,3)==1
                    mod_entry=[mod_entry; temp; squeeze(M(dir_idx.offd_col,:,:))'; squeeze(M(dir_idx.offd_row,:,:))'];
                else
                    mod_entry=[mod_entry; temp; squeeze(M(dir_idx.offd_col,:,:)); squeeze(M(dir_idx.offd_row,:,:))];
                end
                
            else                            % diagonal direction
                mod_entry=squeeze(M(:,dir_idx.diag_row,:));
                mod_entry(dir_idx.diag_row,:)=[];
                temp=squeeze(M(dir_idx.diag_row,:,:));
                if size(M,3)==1
                    mod_entry=[mod_entry; temp'];
                else
                    mod_entry=[mod_entry; temp];
                end
            end
end
