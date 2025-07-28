function [u, v] = KS_onestep_ETD_second(elem, node, dt, uold, vold, u1, v1, fu, fv, M, A, Ac)
    global kappa;
    kappa = 3;
    % 计算中间值 vh
    vh = (vold + v1) / 2;
    
    % 组装 K2 和人工扩散矩阵 D
    K2 = assemblematrix2D(node, elem, vh);
    D = -K2;
    DT = D.';
    D((D - DT) < 0) = D((D - DT) < 0) * 0 + DT((D - DT) < 0);
    D = D - spdiags(diag(D), 0, size(K2, 1), size(K2, 2));
    D(D < 0) = 0;
    D = D - spdiags(sum(D, 2), 0, size(D, 1), size(D, 2));
    
    % 计算 AA
    AA = M \ (K2 + D - A) - kappa * speye(size(K2));
    
    % ETD1 求解 u, v
    u = phipm(dt, AA, [uold, kappa * uold + fu, 1/dt * (kappa * u1 - kappa * uold) + 1/dt * fu], eps, false);
    v = phipm(dt, Ac, [vold, uold + fv, 1/dt * (u1 - uold) + 1/dt * fv], eps, true);
end