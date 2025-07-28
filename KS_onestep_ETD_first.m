function [u1, v1, fu, fv, M, A, Ac] = KS_onestep_ETD_first(elem, node, t, dt, uold, vold, pdeU, pdeV)
    % 组装矩阵
    global alpha kappa;
    alpha = 1; kappa = 3;
    [A,M,~] = assemblematrix(node,elem,1); %质量集中
    N = size(node, 1);
    M = spdiags(M, 0, N, N);
    
    % 组装 K1 和人工扩散矩阵 D1
    K1 = assemblematrix2D(node, elem, vold);
    D1 = -K1;
    DT1 = D1.';
    D1((D1 - DT1) < 0) = D1((D1 - DT1) < 0) * 0 + DT1((D1 - DT1) < 0);
    D1 = D1 - spdiags(diag(D1), 0, size(K1, 1), size(K1, 2));
    D1(D1 < 0) = 0;
    D1 = D1 - spdiags(sum(D1, 2), 0, size(D1, 1), size(D1, 2));
    
    % 计算 AA1 和 Ac
    AA1 = M \ (K1 + D1 - A) - kappa * speye(size(K1));
    Ac = -M \ A - alpha * speye(size(K1));
    
    % 组装右端项 fu, fv
    if isfield(pdeU, 'f') && isfield(pdeV, 'f')
        fu = assembleRHSt3D(node, elem, pdeU, t, []);
        fv = assembleRHSt3D(node, elem, pdeV, t, []);
    else
        fu = 0; fv = 0;
    end
    
    % ETD1 求解 u1, v1
    u1 = phipm(dt, AA1, [uold, fu + kappa * uold], eps, false);
    v1 = phipm(dt, Ac, [vold, fv + uold], eps, true);
end