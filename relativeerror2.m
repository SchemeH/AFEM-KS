function Rerror = relativeerror2(node, elem, error, vh)
% Compute relative error with protection for near-zero solutions

[Dphi, area] = gradbasis(node, elem); 

% Compute numerical solution gradient
dudx = vh(elem(:, 1)) .* squeeze(Dphi(:, 1, 1)) + ...
       vh(elem(:, 2)) .* squeeze(Dphi(:, 1, 2)) + ...
       vh(elem(:, 3)) .* squeeze(Dphi(:, 1, 3));
dudy = vh(elem(:, 1)) .* squeeze(Dphi(:, 2, 1)) + ...
       vh(elem(:, 2)) .* squeeze(Dphi(:, 2, 2)) + ...
       vh(elem(:, 3)) .* squeeze(Dphi(:, 2, 3));

solnDu = [dudx, dudy];
Du_h1 = sum(solnDu.^2, 2);
Du_h1 = sqrt(sum(Du_h1 .* area)); % H1 semi-norm

% Calculate L2 norm of the solution for additional protection
mass = accumarray(elem(:), repmat(area, 3, 1), [size(vh, 1), 1]);
u_L2 = sqrt(sum(vh.^2 .* mass));

% Determine appropriate denominator
if Du_h1 > 1e-8
    denominator = Du_h1; % Use H1 semi-norm when significant
elseif u_L2 > 1e-8
    denominator = u_L2; % Fallback to L2 norm
else
    denominator = 1; % Use absolute error for near-zero solutions
    fprintf('Warning: Solution norm extremely small, using absolute error\n');
end

Rerror = error / denominator; % Relative error
end