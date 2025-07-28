function [eta, Du] = estimaterecovery2(node, elem, u)
%% ESTIMATERECOVERY recovery type error estimator with zero-value suppression
%  
% Modified to suppress error estimation in regions where u is close to zero

[Dlambda, area] = gradbasis(node, elem);
Du = gradu(node, elem, u, Dlambda);
Du = recovery(node, elem, Du, area); % Gradient reconstruction

% Calculate element-wise average of u
u_elem = mean(u(elem), 2); % Average of u at three nodes per element

% Initialize DDu
DDu = zeros(size(elem, 1), 4);

% Only compute second derivatives where u is significant
threshold = 1e-8; % Threshold for considering u as non-zero
idx_nonzero = find(u_elem > threshold);

if ~isempty(idx_nonzero)
    % Compute second derivatives only for non-zero elements
    DDu(idx_nonzero, 1:2) = gradu(node, elem(idx_nonzero, :), Du(idx_nonzero, 1), Dlambda(idx_nonzero, :, :));
    DDu(idx_nonzero, 3:4) = gradu(node, elem(idx_nonzero, :), Du(idx_nonzero, 2), Dlambda(idx_nonzero, :, :));
end

% Calculate error estimator with suppression in zero regions
eta = area .* sum(abs(DDu), 2);

% Apply suppression factor in low-value regions
suppression_factor = min(1, u_elem / (threshold * 10)); % Scale from 0 to 1
eta = eta .* suppression_factor;

% TODO: add diffusion coefficient if needed
end