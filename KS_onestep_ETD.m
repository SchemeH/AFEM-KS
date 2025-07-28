function [u, v, u1, v1] = KS_onestep_ETD(elem, node, t, dt, uold, vold,  pdeU, pdeV)
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
         
        AA1 = M\(K1+D1-A)-kappa*speye(size(K1));
        Ac = -M\A - alpha*speye(size(K1));

        if isfield(pdeU, 'f') && isfield(pdeV, 'f')
             fu  = assembleRHSt2D(node, elem, pdeU, t, []);
             fv  = assembleRHSt2D(node, elem, pdeV, t, []);
        else
             fu = 0; fv = 0;
        end

         % ETD1 求解u
         u1 = phipm(dt,AA1,[uold,fu+kappa*uold],eps,false);
         % u1 = phipm(dt,AA1,[uold,uold],eps,false);
         v1 = phipm(dt,Ac,[vold,uold+fv],eps,true);
         
         % v1 = phipm(dt,Ac,[vold,fv+uold],eps,true);
         vh = vold/2 + v1/2;
          %% 组装左端矩阵
         K2 = assemblematrix2D(node, elem, vh);
         % 组装人工扩散矩阵
         D = -K2;
         DT = D.';
         D((D -(DT))<0) = D((D -(DT))<0)*0  + DT((D -DT)<0);

         D = D - spdiags(diag(D), 0, size(K2, 1), size(K2, 2));
         D(D<0) = 0;
         D = D - spdiags(sum(D, 2), 0, size(D, 1), size(D, 2));
         AA = M\(K2+D-A)-kappa*speye(size(K2));
         

         % ETD1 求解u
         u = phipm(dt,AA,[uold,kappa*uold+fu,1/dt*(kappa*u1 - kappa*uold)+1/dt*fu],eps,false);
         % u = phipm(dt,AA,[uold,uold,1/dt*(u1 - uold)],eps,false);
         % v = phipm(dt,Ac,[vold,uold+fv,1/dt*(u1 - uold)+1/dt*fv],eps,true);
         v = phipm(dt,Ac,[vold,uold+fv,1/dt*(u1 - uold)+1/dt*fv],eps,true);
end