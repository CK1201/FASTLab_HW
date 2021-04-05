function Ct = getCt(n_seg, n_order)
    %#####################################################
    % STEP 2.1: finish the expression of Ct
    d_size=(n_order+1)/2;
    d_all_size = 2*n_seg*d_size;
    dFP_size=(n_seg+1)*d_size;
    Ct=zeros(d_all_size,dFP_size);
    % start end
    for i=1:d_size
        Ct(i,i)=1;
        Ct(d_all_size-d_size+i,d_size+n_seg-1+i)=1;
    end
    % waypoint
    for i=1:n_seg-1
        Ct(d_size+2*(i-1)*d_size+1,d_size+i)=1;% p_end
        Ct(d_size+2*(i-1)*d_size+1+d_size,d_size+i)=1;% p_start

        Ct(d_size+2*(i-1)*d_size+2,2*d_size+n_seg+(i-1)*(d_size-1)+0)=1;% v_end
        Ct(d_size+2*(i-1)*d_size+3,2*d_size+n_seg+(i-1)*(d_size-1)+1)=1;% a_end
        Ct(d_size+2*(i-1)*d_size+4,2*d_size+n_seg+(i-1)*(d_size-1)+2)=1;% j_end

        Ct(d_size+2*(i-1)*d_size+2+d_size,2*d_size+n_seg+(i-1)*(d_size-1)+0)=1;% v_start
        Ct(d_size+2*(i-1)*d_size+3+d_size,2*d_size+n_seg+(i-1)*(d_size-1)+1)=1;% a_start
        Ct(d_size+2*(i-1)*d_size+4+d_size,2*d_size+n_seg+(i-1)*(d_size-1)+2)=1;% j_start
    end
    
end