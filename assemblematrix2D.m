function K = assemblematrix2D(node, elem, vh)
%% vh是定义在二维网格上的数值解，(N,1)数组，N为节点数
N = size(node, 1);
NT = size(elem, 1);

% 计算基函数梯度和面积（二维用面积代替体积）
[Dphi, area] = gradbasis(node,elem); % 改为二维梯度计算函数（NT,2,3）
dudx = vh(elem(:,1)).*Dphi(:,1,1) + vh(elem(:,2)).*Dphi(:,1,2) + vh(elem(:,3)).*Dphi(:,1,3);
dudy = vh(elem(:,1)).*Dphi(:,2,1) + vh(elem(:,2)).*Dphi(:,2,2) + vh(elem(:,3)).*Dphi(:,2,3); 
solnDu = [dudx, dudy]; % 改为二维梯度（NT,2）

% 获取二维积分点和权重（例如4阶积分）
[lambda, weight] = quadpts(4); % 改为二维积分函数
t1 = sum( (weight.') .* lambda, 1); % (1,3) 每个基函数的积分权重
basis_integrals = area.*t1; % (NT,3) 每个单元基函数的积分

%% 组装K矩阵：K(i,j) = ∫φ_j(∇vh·∇φ_i)dx
K = sparse(N, N);
for i = 1:3 % 三角形三个顶点
    for j = 1:3
        ii = double(elem(:,i)); % 当前单元第i个节点的全局索引
        jj = double(elem(:,j));
        % 二维梯度点乘：∇vh·∇φ_i = dudx*Dphi_x + dudy*Dphi_y
        Kij = (solnDu(:,1).*Dphi(:,1,i) + solnDu(:,2).*Dphi(:,2,i)) .* basis_integrals(:,j);
        K = K + sparse(ii, jj, Kij, N, N);
    end
end
end