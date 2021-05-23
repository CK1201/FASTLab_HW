function [Aieq, bieq] = getAbieq(n_seg, n_order, corridor_range, ts, v_max, a_max)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % STEP 3.2.1 p constraint
    Aieq_p_k = eye(n_all_poly);
    Aieq_p = [-Aieq_p_k; Aieq_p_k];
    bieq_p = zeros(2*n_all_poly,1);
    for i = 1:n_seg
        bieq_p((i-1)*(n_order+1)+1:i*(n_order+1)) = -corridor_range(i,1);
        bieq_p(n_all_poly + (i-1)*(n_order+1)+1:n_all_poly + i*(n_order+1)) = corridor_range(i,2);
    end
    
    % connection
    Aieq_p_con = zeros(2 * n_seg, n_all_poly);
    bieq_p_con = zeros(2 * n_seg,1);
    for i = 1:n_seg - 1
        Aieq_p_con(i, i * (n_order+1)) = -1;
        Aieq_p_con(i + n_seg, i * (n_order+1)) = 1;
        corridor_range_temp = [corridor_range(i,:),corridor_range(i+1,:)];
        corridor_range_temp = sort(corridor_range_temp);
%         disp(corridor_range_temp)
        bieq_p_con(i) = -corridor_range_temp(2);
        bieq_p_con(i + n_seg) = corridor_range_temp(3);
    end
    Aieq_p = [Aieq_p; Aieq_p_con];
    bieq_p = [bieq_p; bieq_p_con];
    %#####################################################
    % STEP 3.2.2 v constraint 
    Aieq_v = [];
    bieq_v = [];
    Aieq_v = zeros(2*(n_all_poly-1),n_all_poly);
    bieq_v = zeros(2*(n_all_poly-1),1);
    for i = 1:n_all_poly-1
        Aieq_v(i, i+1) = -n_order;
        Aieq_v(i, i) = n_order;
        bieq_v(i) = v_max;
        
        Aieq_v(i+n_all_poly-1, i+1) = n_order;
        Aieq_v(i+n_all_poly-1, i) = -n_order;
        bieq_v(i+n_all_poly-1) = v_max;
    end

    %#####################################################
    % STEP 3.2.3 a constraint   
    Aieq_a = [];
    bieq_a = [];
    Aieq_v = zeros(2*(n_all_poly-2),n_all_poly);
    bieq_v = zeros(2*(n_all_poly-2),1);
    for i = 1:n_all_poly-2
        Aieq_v(i+n_all_poly-1, i+2) = -n_order*(n_order-1);
        Aieq_v(i+n_all_poly-1, i+1) = 2*n_order*(n_order-1);
        Aieq_v(i+n_all_poly-1, i) = -n_order*(n_order-1);
        bieq_v(i) = a_max;
        
        Aieq_v(i+n_all_poly-1, i+2) = n_order*(n_order-1);
        Aieq_v(i+n_all_poly-1, i+1) = -2*n_order*(n_order-1);
        Aieq_v(i+n_all_poly-1, i) = n_order*(n_order-1);
        bieq_v(i+n_all_poly-1) = a_max;
    end
    
    %#####################################################
    % combine all components to form Aieq and bieq   
    Aieq = [Aieq_p; Aieq_v; Aieq_a];
    bieq = [bieq_p; bieq_v; bieq_a];
%     Aieq = Aieq_p;
%     bieq = bieq_p;
end