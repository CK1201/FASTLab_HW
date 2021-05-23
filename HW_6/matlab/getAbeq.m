function [Aeq, beq] = getAbeq(n_seg, n_order, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    M=[];
    M_k = getM(n_order);
    for k = 1:n_seg
        M = blkdiag(M, M_k);
    end
    %#####################################################
    % STEP 2.1 p,v,a constraint in start 
    Aeq_start = zeros(4, n_all_poly);
    beq_start = zeros(4, 1);
    % STEP 2.1: write expression of Aeq_start and beq_start
    Aeq_start(1,1)=1;% p
    Aeq_start(2,2)=1;% v
    Aeq_start(3,3)=2;% a
    Aeq_start(4,4)=6;% j
    for i=1:length(start_cond)
        beq_start(i)=start_cond(i);
    end

    %#####################################################
    % STEP 2.2 p,v,a constraint in end
    Aeq_end = zeros(4, n_all_poly);
    beq_end = zeros(4, 1);
    % STEP 2.2: write expression of Aeq_end and beq_end
    % p
    for j=0:n_order
            Aeq_end(1,n_all_poly-n_order+j)=ts(end)^j;
    end
    % v
    for j=1:n_order
        Aeq_end(2,n_all_poly-n_order+j)=j*ts(end)^(j-1);
    end
    % a
    for j=2:n_order
        Aeq_end(3,n_all_poly-n_order+j)=j*(j-1)*ts(end)^(j-2);
    end
    % j
    for j=3:n_order
        Aeq_end(4,n_all_poly-n_order+j)=j*(j-1)*(j-2)*ts(end)^(j-3);
    end

    for i=1:length(end_cond)
        beq_end(i)=end_cond(i);
    end
    
    %#####################################################
    % position constrain in all middle waypoints
%     Aeq_wp = zeros(n_seg-1, n_all_poly);
%     beq_wp = zeros(n_seg-1, 1);
%     % STEP 2.3: write expression of Aeq_wp and beq_wp
%     % end of seg
%     for i=1:n_seg-1
%         for j=0:n_order
%             Aeq_wp(i,i*(n_order+1)-n_order+j)=ts(i)^j;
%         end
%         beq_wp(i)=waypoints(i+1);
%     end

    %#####################################################
    % STEP 2.3 position continuity constrain between 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);
    for i=1:n_seg-1
        for j=0:n_order
            Aeq_con_p(i,i*(n_order+1)-n_order+j)=ts(i)^j;
        end
        Aeq_con_p(i,i*(n_order+1)+1)=-1;
    end

    %#####################################################
    % STEP 2.4 velocity continuity constrain between 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);
    for i=1:n_seg-1
        for j=1:n_order
            Aeq_con_v(i,i*(n_order+1)-n_order+j)=j*ts(i)^(j-1);
        end
        Aeq_con_v(i,i*(n_order+1)+2)=-1;
    end

    %#####################################################
    % STEP 2.5 acceleration continuity constrain between 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);
    for i=1:n_seg-1
        for j=2:n_order
            Aeq_con_a(i,i*(n_order+1)-n_order+j)=j*(j-1)*ts(i)^(j-2);
        end
        Aeq_con_a(i,i*(n_order+1)+3)=-2;
    end

    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a];
    beq_con = [beq_con_p; beq_con_v; beq_con_a];
    Aeq = [Aeq_start; Aeq_end; Aeq_con];
    Aeq = Aeq * M;
    beq = [beq_start; beq_end; beq_con];
end