function [u, v] =KS_onestep_CN(elem, node, t, dt, uold, vold,  pdeU, pdeV)
         [A,M,~] = assemblematrix(node,elem,1); %��������
         N = size(node, 1);
         M =  spdiags(M,0,N,N);
        %% ��װ��˾���
         K = assemblematrix2D(node, elem, vold);
         % ��װ�˹���ɢ����
         D = -K;
         DT = D.';
         D((D -(DT))<0) = D((D -(DT))<0)*0  + DT((D -DT)<0);

         D = D - spdiags(diag(D), 0, size(K, 1), size(K, 2));
         D(D<0) = 0;
         D = D - spdiags(sum(D, 2), 0, size(D, 1), size(D, 2));
         AA = M + dt/2*(A - K - D);
         
        % ��װ�Ҷ���
        if isfield(pdeU, 'f') && isfield(pdeV, 'f')
             fu  = assembleRHSt2D(node, elem, pdeU, t, []);
             fv  = assembleRHSt2D(node, elem, pdeV, t, []);
             size(M)
             size(uold)
             size(fu)
             bu = (M - dt/2*(A - K - D))*uold + dt*fu;
        else
             bu = (M - dt/2*(A - K - D))*uold;
        end
         
         % ���u
         %u = AA\bu;
         if N > 50000
            u = amg(AA,bu);
         else
             u = AA\bu;
         end
         %% ���v
         if isfield(pdeU, 'f') && isfield(pdeV, 'f')
             F = dt/2*M*(u + uold);
             bv = (M - dt/2*A - dt/2*M)*vold + F + dt * fv;
         else
             F = dt/2*M*(u + uold);
             bv = (M - dt/2*A - dt/2*M)*vold + F;
         end
         AA = M + dt/2*A + dt/2*M;
         %v = AA\bv;
         if N > 50000
            v = amg(AA,bv);
         else
            v = AA\bu;
         end
end