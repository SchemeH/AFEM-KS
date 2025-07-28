function pde = KSdataconvergence(Te)
    % 定义 pdeU 的结构体
    pde = struct('f',@f,'exactu',@exactu,'g_Du',@g_Du,'Du',@Du,'u',@u,'du',@du,'w',@w,'exactv',@exactv,'g_Dv',@g_Dv,'Dv',@Dv,'v',@v,'dv',@dv);
    % pdeU 的右端项函数
    function s = f(p, t) 
        % 右端项
        x = p(:, 1); y = p(:,2);
        s = (3 - 4*(x.^2 + y.^2)).*exp(-(x.^2 + y.^2 + t)) + (-4 + 8*(x.^2 + y.^2)).*exp(-2*(x.^2+ y.^2 + t));
    end

    % pdeU 的 Dirichlet 边界条件函数
    function s = g_Du(p, t) 
        s = exactu(p, t);
    end

    % pdeU 的解析解函数
    function z = exactu(p, t)
        % 解析解
        x = p(:,1); y = p(:,2);
        z = exp(-(x.^2 + y.^2 + t));
    end

    % pdeU 的 Du 函数
    function z = Du(p, t)    
        % 解析解
        x = p(:,1); y = p(:,2); 
        z(:,1) = exp(-(x.^2 + y.^2 + t)).*(-2*x);
        z(:,2) = exp(-(x.^2 + y.^2 + t)).*(-2*y);
    end

    % pdeU 的 u 函数
    function s = u(p)
        s = pde.exactu(p, Te);
    end

    % pdeU 的 du 函数
    function s = du(p)
        s = pde.Du(p, Te);
    end



    % pdeV 的右端项函数
    function s = w(p, t) 
        % 右端项
        x = p(:,1); y = p(:,2);
        s = (3 - 4*(x.^2 + y.^2)).*exp(-(x.^2 + y.^2 + t));
    end

    % pdeV 的 Dirichlet 边界条件函数
    function s = g_Dv(p, t) 
        s = exactv(p, t);
    end

    % pdeV 的解析解函数
    function z = exactv(p, t)
        % 解析解
        x = p(:,1); y = p(:,2);
        z = exp(-(x.^2 + y.^2 + t));
    end

    % pdeV 的 Du 函数
    function z = Dv(p, t)    
        % 解析解
        x = p(:,1); y = p(:,2); 
        z(:,1) = exp(-(x.^2 + y.^2 + t)).*(-2*x);
        z(:,2) = exp(-(x.^2 + y.^2 + t)).*(-2*y);
    end

    % pdeV 的 u 函数
    function s = v(p)
        s = pde.exactv(p, Te);
    end

    % pdeV 的 du 函数
    function s = dv(p)
        s = pde.Dv(p, Te);
    end
end