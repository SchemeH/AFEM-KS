function markedElem = mark2(elem, eta, theta, method)
% MARK2 Enhanced marking function with combined refinement/coarsening
%
% markedElem = mark2(elem,eta,theta) mark elements using Dorfler strategy.
%
% New feature: 'HYBRID' method for combined refinement and coarsening
%   - Marks high-gradient areas for refinement (thickening)
%   - Marks all other areas for coarsening
%
% Additional methods:
%   'REFINE' - Marks high-gradient areas (eta > theta*max(eta))

NT = size(elem,1); 
isMark = false(NT,1);
if ~exist('method','var'), method = 'L2'; end  % default marking

switch upper(method)
    case 'MAX'
        isMark(eta > theta*max(eta)) = 1;
        
    case 'REFINE'  % 新增：专门标记高梯度区域
        isMark(eta > theta*max(eta)) = 1;
        fprintf('[REFINE] 标记比例: %.1f%%\n', 100*nnz(isMark)/NT);
        
    case 'COARSEN'
        % 原始粗化策略
        isMark(eta < theta*max(eta)) = 1;
        
    case 'HYBRID'  % 新增：混合策略（核心修改）
        % ===== 第一步：标记高梯度区域（用于细化/加粗）=====
        refineThreshold = theta * max(eta);
        toRefine = eta > refineThreshold;
        
        % ===== 第二步：标记所有其他区域（用于粗化）=====
        toCoarsen = ~toRefine;  % 所有非高梯度区域
        
        % ===== 组合标记 =====
        isMark = toRefine | toCoarsen;  % 所有元素都被标记
        
        % 诊断信息
        fprintf('[HYBRID] 细化区域: %.1f%%, 粗化区域: %.1f%%\n', ...
                100*nnz(toRefine)/NT, 100*nnz(toCoarsen)/NT);
        
    case 'L2'
        [sortedEta,idx] = sort(eta.^2,'descend'); 
        x = cumsum(sortedEta);
        isMark(idx(x < theta* x(NT))) = 1;
        isMark(idx(1)) = 1;
end

markedElem = uint32(find(isMark));
end