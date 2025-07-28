function pde = KSdataV(T)
pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du,'u',@u,'du',@du);
    
    function s = f(p, t) % load data (right hand side function)
    % ”“∂ÀœÓ
    x = p(:,1); y = p(:,2);
    s = (3 - 4*(x.^2 + y.^2)).*exp(-(x.^2 + y.^2 + t));
    end

    function s = g_D(p, t) % Dirichlet boundary condition
    s = exactu(p, t);
    end

    function z=exactu(p, t)
    %Ω‚ŒˆΩ‚
        x = p(:,1); y = p(:,2);
        z = exp(-(x.^2 + y.^2 + t));
    end

    function z = Du(p, t)    
    % Ëß£ÊûêËß?
        x = p(:,1); y = p(:,2); 
        z(:,1) = exp(-(x.^2 + y.^2 + t)).*(-2*x);
        z(:,2) = exp(-(x.^2 + y.^2 + t)).*(-2*y);
    end

    function s = u(p)
    s = pde.exactu(p, T);
    end

    function s = du(p)
    s = pde.Du(p, T);
    end

end