function M = getM(n_seg, n_order, ts)
    M = [];
    for k = 1:n_seg
%         M_k = [];
        %#####################################################
        % STEP 1.1: calculate M_k of the k-th segment 
        M_k=zeros(n_order+1,n_order+1);
        M_k(1,1)=1;
        M_k(2,2)=1;
        M_k(3,3)=2;
        M_k(4,4)=6;
        % p
        for j=0:n_order
            M_k(5,j+1)=ts(end)^j;
        end
        % v
        for j=1:n_order
            M_k(6,j+1)=j*ts(end)^(j-1);
        end
        % a
        for j=2:n_order
            M_k(7,j+1)=j*(j-1)*ts(end)^(j-2);
        end
        % j
        for j=3:n_order
            M_k(8,j+1)=j*(j-1)*(j-2)*ts(end)^(j-3);
        end
        
        M = blkdiag(M, M_k);
    end
end