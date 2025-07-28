function pde = heatafemdata(T)
%% SQUAREAFEMDATA Data of an example of AFEM in a square
%
%  u = \beta(t) exp(-[(x - t + 0.5)^2 + (y - t + 0.5)^2]/0.04 with,
%  \beta(t) = 0.1(1 - exp(-10^2(t - 0.5)^2)).
%  
pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du,'u',@u,'du',@du);
    
    function s = f(p, t) % load data (right hand side function)
    % ”“∂ÀœÓ
    x = p(:,1); y = p(:,2);
    beta = 0.1*(1 - exp(-100*(t - 0.5).^2));
    u1 = exp(   -(  (x - t + 0.5 ).^2  +  (y - t + 0.5).^2  ) /0.04  );
    betaDt = 20*exp(-100*(t - 0.5).^2).*(t - 0.5);
    u1Dt = u1.*(  x +y - 2*t  + 1 )/0.02;
    u_t = betaDt.*u1 + u1Dt.*beta;
    
    u_x = -(x - t + 0.5)/0.02.*beta.*u1;
    u_y = -(y - t + 0.5)/0.02.*beta.*u1;
    u_xx = -1/0.02*beta.*u1 - 1/0.02*(x - t + 0.5).*u_x;
    u_yy = -1/0.02*beta.*u1 - 1/0.02*(y - t + 0.5).*u_y;
    s = u_t - u_xx - u_yy;
    end

    function s = g_D(p, t) % Dirichlet boundary condition
    s = exactu(p, t);
    end

    function z=exactu(p, t)
    %Ω‚ŒˆΩ‚
        x = p(:,1); y = p(:,2);
        beta = 0.1*(1 - exp(-100*(t - 0.5).^2));
        z = beta.*exp(   -(  (x - t + 0.5 ).^2  +  (y - t + 0.5).^2) /0.04  );
    end

    function z = Du(p, t)    
    % Ëß£ÊûêËß?
        x = p(:,1); y = p(:,2); 
        beta = 0.1*(1 - exp(-100*(t - 0.5)^2));
        u1 = exp(   -(  (x - t + 0.5 ).^2  +  (y - t + 0.5).^2  ) /0.04  );  
        u_x = -(x - t + 0.5)/0.02.*(beta.*u1);
        u_y = -(y - t + 0.5)/0.02.*(beta.*u1);
        z(:,1) = u_x;
        z(:,2) = u_y;
    end

    function s = u(p)
    s = pde.exactu(p, T);
    end

    function s = du(p)
    s = pde.Du(p, T);
    end

end