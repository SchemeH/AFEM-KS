function [u, v] = KS_onestep_ETD1(elem, node, t, dt, uold, vold,  pdeU, pdeV)
        global kappa alpha;
        kappa = 0; alpha = 1;
        [A,M,~] = assemblematrix(node,elem,1); %质量集中
        N = size(node, 1);
        M =  spdiags(M,0,N,N);
        %% 组装左端矩阵
        K1 = assemblematrix2D(node, elem, vold);
        % 组装人工扩散矩阵
        D1 = -K1;
        DT1 = D1.';
        D1((D1 -(DT1))<0) = D1((D1 -(DT1))<0)*0  + DT1((D1 -DT1)<0);

        D1 = D1 - spdiags(diag(D1), 0, size(K1, 1), size(K1, 2));
        D1(D1<0) = 0;
        D1 = D1 - spdiags(sum(D1, 2), 0, size(D1, 1), size(D1, 2));
         
        AA = M\(K1+D1-A-kappa*speye(size(K1)));
        Ac = -M\A - alpha*speye(size(K1));

        if isfield(pdeU, 'f') && isfield(pdeV, 'f')
             fu  = assembleRHSt2D(node, elem, pdeU, t, []);
             fv  = assembleRHSt2D(node, elem, pdeV, t, []);
        else
             fu = 0; fv = 0;
        end

         % ETD1 求解u
         u = phipm(dt,AA,[uold,fu+kappa*uold],eps,false);
         v = phipm(dt,Ac,[vold,fv+uold],eps,true);
end