function L2error = timeEstimateL2(node, elem, uh, uhold)
    delta_u = uh - uhold;
    [lambda, weight] = quadpts(3); % 二维积分点（可能为三角形7点积分）
    
    % 调整节点索引为三角形单元（3个顶点）
    elem2uh = [delta_u(elem(:,1)), delta_u(elem(:,2)), delta_u(elem(:,3))]; % (NT,3)
    
    % 计算每个积分点的函数值（注意维度匹配）
    uh_value = lambda*elem2uh'; % (NQ,NT)
    
    % 计算平方积分（二维权重要素处理）
    eta = sum(uh_value.^2 .* weight(:), 1); % (1,NT)
    
    % 计算三角形面积（替换三维体积计算）
    d12 = node(elem(:,2),:) - node(elem(:,1),:);
    d13 = node(elem(:,3),:) - node(elem(:,1),:);
    area = abs(d12(:,1).*d13(:,2) - d12(:,2).*d13(:,1))/2; % 叉积计算面积
    
    % 最终L2误差
    L2error = sum(area(:) .* eta(:));
end