function Rerror = relativeerror(node, elem, error, vh)
% 计算二维相对误差
    [Dphi, area] = gradbasis(node, elem);  % 基函数梯度 Dphi (NT, 2, 3)
    
    % 计算数值解的梯度场
    dudx = vh(elem(:,1)).*Dphi(:,1,1) + vh(elem(:,2)).*Dphi(:,1,2)...
          + vh(elem(:,3)).*Dphi(:,1,3);
    dudy = vh(elem(:,1)).*Dphi(:,2,1) + vh(elem(:,2)).*Dphi(:,2,2)...
          + vh(elem(:,3)).*Dphi(:,2,3);
    
    solnDu = [dudx, dudy];       % 数值解梯度矩阵 (NT, 2)
    Du_h1 = sum(solnDu.^2, 2);   % 每个单元梯度模方
    Du_h1 = sum(Du_h1.*area);    % 全域H1半模
    Du_h1 = sqrt(Du_h1);         % H1范数
    Rerror = error/Du_h1;        % 相对误差
end