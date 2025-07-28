function pde = KSdataU(T)


pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du,'u',@u,'du',@du);
    
    function s = f(p, t) % load data (right hand side function)
    % ”“∂ÀœÓ
    x = p(:, 1); y = p(:,2);
    % s = (3 - 4*(x.^2 + y.^2)).*exp(-(x.^2 + y.^2 + t)) + (-4 + 8*(x.^2 + y.^2)).*exp(-2*(x.^2+ y.^2 + t));
    s = -2*exp(-2*t).*(cos(x).^2 + cos(x).*cos(y) + cos(y).^2 - 1);
    end

    function s = g_D(p, t) % Dirichlet boundary condition
    s = exactu(p, t);
    end

    function z=exactu(p, t)
    %Ω‚ŒˆΩ‚
        x = p(:,1); y = p(:,2);
        % z = exp(-(x.^2 + y.^2 + t));
        % z = exp(t).*((x.^2 - x).^2).*(y.^2 - y).^2;
        z = exp(-t).*(cos(x) + cos(y));
    end

    function z = Du(p, t)    
    % Ëß£ÊûêËß?
        x = p(:,1); y = p(:,2); 
        % z(:,1) = exp(-(x.^2 + y.^2 + t)).*(-2*x);
        % z(:,2) = exp(-(x.^2 + y.^2 + t)).*(-2*y);
        z(:,1) = 2*exp(t).*((y.^2 -y).^2).*(x.^2 -x).*(2*x-1);
        z(:,2) = 2*exp(t).*((x.^2 -x).^2).*(y.^2 -y).*(2*y-1);
    end

    function s = u(p)
    s = pde.exactu(p, T);
    end

    function s = du(p)
    s = pde.Du(p, T);
    end

end