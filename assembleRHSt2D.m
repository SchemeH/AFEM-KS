function b = assembleRHSt2D(node, elem, pde, t, option)
N = size(node,1); 
NT = size(elem,1);
Ndof = N;
b = zeros(Ndof,1);
[~,area] = gradbasis(node,elem); % 确保使用二维梯度与面积计算
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % 默认积分阶数
end
if ~isempty(pde.f) && ~isreal(pde.f)  % f是函数句柄
    [lambda,weight] = quadpts(option.fquadorder); % 二维积分点
    phi = lambda;                 % 线性基函数
    nQuad = size(lambda,1);
    bt = zeros(NT, 3);            % 每个三角形3个节点
    for p = 1:nQuad
        % 二维积分点坐标计算（使用三角形单元）
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy, t);       % 计算源项
        for i = 1:3
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        end
    end
    bt = bt.*repmat(area,1,3);    % 乘以三角形面积
    b = accumarray(elem(:), bt(:), [Ndof 1]); % 组装全局向量
end
end