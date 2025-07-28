function [u, v] = KS_onestep(elem, node, dt, uold, vold)
         [A,M,~] = assemblematrix(node,elem,1); %质量集中
         N = size(node, 1);
         M =  spdiags(M,0,N,N);
       
         bu = M*uold;
         %% 组装左端矩阵
         K = assemblematrix2D(node, elem, vold);
         % 组装人工扩散矩阵
         D = -K;
         DT = D.';
         D((D -(DT))<0) = D((D -(DT))<0)*0  + DT((D -DT)<0);

         D = D - spdiags(diag(D), 0, size(K, 1), size(K, 2));
         D(D<0) = 0;
         D = D - spdiags(sum(D, 2), 0, size(D, 1), size(D, 2));
         
        
         
         AA = M - dt*( -A + K + D );
         
         % 求解u
         u = AA\bu;
         %% 求解v
         bv = M * vold;
         AA = M + dt * (A + M);
         v = AA\bv;
end