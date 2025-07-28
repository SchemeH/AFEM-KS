function [eta, Du, eta_r] = relativeerror3(node,elem,u,normType)
%% ESTIMATERECOVERY recovery type error estimator with selectable norm type.
%  
% [eta, Du, eta_r] = estimaterecovery(node,elem,u) computes:
%   eta: error estimator 
%   Du: recovered gradient
%   eta_r: relative error in selected norm
%
% [eta, Du, eta_r] = estimaterecovery(node,elem,u,normType) specifies the norm type:
%   'L2' - L² norm (default)
%   'max' - maximum norm
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

if nargin < 4
    normType = 'L2'; % 默认使用L²范数
end

% 预计算梯度基函数和面积
[Dlambda,area] = gradbasis(node,elem);
% 保存原始梯度 (分片常数)
original_grad = gradu(node,elem,u,Dlambda); 
% 恢复梯度 (P1 连续)
recovered_grad = recovery(node,elem,original_grad,area); 
Du = recovered_grad; % 输出恢复梯度

% 初始化误差指示器和相对误差
eta = zeros(size(elem,1),1);
eta_r = 0;

% 计算每个单元的L2误差和恢复梯度范数
num_squares = zeros(size(elem,1),1); % |R_h - ∇u_h|² per element
den_squares = 0; % |R_h|² globally

for e = 1:size(elem,1)
    areaK = area(e);
    nodes = elem(e,:);
    gx = recovered_grad(nodes,1);
    gy = recovered_grad(nodes,2);
    cKx = original_grad(e,1);
    cKy = original_grad(e,2);
    
    % 计算单元上恢复梯度与原始梯度的L2差方
    int_diffx2 = areaK * sum( ...
        (gx - cKx).^2 + ...
        (gx - cKx).*circshift(gx - cKx,1) ) / 6;
    int_diffy2 = areaK * sum( ...
        (gy - cKy).^2 + ...
        (gy - cKy).*circshift(gy - cKy,1) ) / 6;
    num_squares(e) = int_diffx2 + int_diffy2;
    
    % 计算单元上恢复梯度的L2范数
    if nargout >= 3 && strcmpi(normType,'l2')
        int_gx2 = areaK * sum(gx.^2 + gx.*circshift(gx,1)) / 6;
        int_gy2 = areaK * sum(gy.^2 + gy.*circshift(gy,1)) / 6;
        den_squares = den_squares + int_gx2 + int_gy2;
    end
end

% 误差指示器: 单元上|R_h - ∇u_h|的L2范数
eta = sqrt(num_squares);

% 计算相对误差
if nargout >= 3
    switch lower(normType)
        case 'l2'
            num = sum(num_squares);
            if den_squares < eps
                eta_r = 0;
            else
                eta_r = sqrt(num) / sqrt(den_squares);
            end
            
        case 'max'
            max_diff = 0;   % ||R_h - ∇u_h||_∞
            recovered_norm = sqrt(recovered_grad(:,1).^2 + recovered_grad(:,2).^2);
            max_recovered = max(recovered_norm);
            
            % 仅在单元顶点处计算差值
            for e = 1:size(elem,1)
                nodes = elem(e,:);
                cKx = original_grad(e,1);
                cKy = original_grad(e,2);
                for i = 1:3
                    node_i = nodes(i);
                    gx = recovered_grad(node_i,1);
                    gy = recovered_grad(node_i,2);
                    diff_norm = sqrt((gx - cKx)^2 + (gy - cKy)^2);
                    if diff_norm > max_diff
                        max_diff = diff_norm;
                    end
                end
            end
            
            if max_recovered < eps
                eta_r = 0;
            else
                eta_r = max_diff / max_recovered;
            end
            
        otherwise
            error('Unsupported norm type. Use ''L2'' or ''max''.');
    end
end
end

% function [eta, Du, eta_r] = relativeerror3(node,elem,u,normType)
% %% ESTIMATERECOVERY recovery type error estimator with selectable norm type.
% %  
% % [eta, Du, eta_r] = estimaterecovery(node,elem,u) computes:
% %   eta: error estimator 
% %   Du: recovered gradient
% %   eta_r: relative error in selected norm
% %
% % [eta, Du, eta_r] = estimaterecovery(node,elem,u,normType) specifies the norm type:
% %   'L2' - L² norm (default)
% %   'max' - maximum norm
% %
% % Copyright (C)  Long Chen. See COPYRIGHT.txt for details.
% 
% % 设置默认范数类型
% if nargin < 4
%     normType = 'L2'; % 默认使用L²范数
% end
% 
% % 预计算节点-单元关联关系（用于最大范数计算）
% node2elem = cell(size(node,1),1);
% for i = 1:size(elem,1)
%     for j = 1:3
%         node2elem{elem(i,j)}(end+1) = i;
%     end
% end
% 
% [Dlambda,area] = gradbasis(node,elem);
% % 保存原始梯度 (分片常数)
% original_grad = gradu(node,elem,u,Dlambda); 
% % 恢复梯度 (P1 连续)
% recovered_grad = recovery(node,elem,original_grad,area); 
% Du = recovered_grad; % 输出恢复梯度
% 
% % 计算原误差估计器
% DDu(:,1:2) = gradu(node,elem,recovered_grad(:,1),Dlambda);
% DDu(:,3:4) = gradu(node,elem,recovered_grad(:,2),Dlambda);
% eta = area.*sum(abs(DDu),2);
% 
% % 计算相对误差 eta_r
% if nargout >= 3
%     switch lower(normType)
%         case 'l2'
%             % L²范数计算
%             num = 0; % 分子累加器 ||R_h - ∇u_h||^2
%             den = 0; % 分母累加器 ||R_h||^2
% 
%             % 计算每个单元的贡献
%             for e = 1:size(elem,1)
%                 areaK = area(e);
%                 nodes = elem(e,:);
% 
%                 % 获取梯度值
%                 gx = recovered_grad(nodes,1);
%                 gy = recovered_grad(nodes,2);
%                 cKx = original_grad(e,1);
%                 cKy = original_grad(e,2);
% 
%                 % 计算恢复梯度的L2范数平方
%                 int_gx2 = areaK * (gx(1)^2 + gx(2)^2 + gx(3)^2 + ...
%                             gx(1)*gx(2) + gx(1)*gx(3) + gx(2)*gx(3)) / 6;
%                 int_gy2 = areaK * (gy(1)^2 + gy(2)^2 + gy(3)^2 + ...
%                             gy(1)*gy(2) + gy(1)*gy(3) + gy(2)*gy(3)) / 6;
%                 den = den + int_gx2 + int_gy2;
% 
%                 % 计算差值的L2范数平方
%                 int_diffx2 = areaK * ( ...
%                     (gx(1)-cKx)^2 + (gx(2)-cKx)^2 + (gx(3)-cKx)^2 + ...
%                     (gx(1)-cKx)*(gx(2)-cKx) + ...
%                     (gx(1)-cKx)*(gx(3)-cKx) + ...
%                     (gx(2)-cKx)*(gx(3)-cKx) ) / 6;
% 
%                 int_diffy2 = areaK * ( ...
%                     (gy(1)-cKy)^2 + (gy(2)-cKy)^2 + (gy(3)-cKy)^2 + ...
%                     (gy(1)-cKy)*(gy(2)-cKy) + ...
%                     (gy(1)-cKy)*(gy(3)-cKy) + ...
%                     (gy(2)-cKy)*(gy(3)-cKy) ) / 6;
% 
%                 num = num + int_diffx2 + int_diffy2;
%             end
% 
%             % 处理分母为零的情况
%             if den < eps
%                 eta_r = 0;
%             else
%                 eta_r = sqrt(num) / sqrt(den);
%             end
% 
%         case 'max'
%             % 最大范数计算
%             max_diff = 0;   % ||R_h - ∇u_h||_∞
%             max_recovered = 0; % ||R_h||_∞
% 
%             % 计算恢复梯度的最大范数
%             recovered_norm = sqrt(recovered_grad(:,1).^2 + recovered_grad(:,2).^2);
%             max_recovered = max(recovered_norm);
% 
%             % 计算差值的最大范数
%             num_nodes = size(node,1);
%             for i = 1:num_nodes
%                 containing_elems = node2elem{i};
% 
%                 for e_idx = 1:length(containing_elems)
%                     e = containing_elems(e_idx);
% 
%                     % 单元上的原始梯度
%                     cKx = original_grad(e,1);
%                     cKy = original_grad(e,2);
% 
%                     % 节点处的恢复梯度
%                     gx = recovered_grad(i,1);
%                     gy = recovered_grad(i,2);
% 
%                     % 计算差值的范数
%                     diff_norm = sqrt((gx - cKx)^2 + (gy - cKy)^2);
% 
%                     % 更新最大值
%                     if diff_norm > max_diff
%                         max_diff = diff_norm;
%                     end
%                 end
%             end
% 
%             % 处理分母为零的情况
%             if max_recovered < eps
%                 eta_r = 0;
%             else
%                 eta_r = max_diff / max_recovered;
%             end
% 
%         otherwise
%             error('Unsupported norm type. Use ''L2'' or ''max''.');
%     end
% end
% end