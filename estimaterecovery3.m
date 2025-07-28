function [eta, Du, eta_r] = estimaterecovery3(node,elem,u)
%% ESTIMATERECOVERY recovery type error estimator.
%  
% eta = estimaterecovery(node,elem,u) computes an error estimator eta by
% recovery second derivative of the gradient of a finite element function u. 
%
% [eta, Du] = estimaterecovery(node,elem,u) also returns the recovered
% derivative Du which is in P1 finite element space.
%
% [eta, Du, eta_r] = estimaterecovery(node,elem,u) additionally returns the
% relative error estimator eta_r = ||R_h(u_h^n) - ∇u_h^n|| / ||R_h(u_h^n)||.
%
% By interpolation error estimate $|u-uI|_{1,2}\leq C|u|_{2,1}$. Therefore
% we recovery an approximation of second derivatives of u and compute L1
% norm. We use the weighted averaging recovery scheme with area weight.
%
% See also recovery, estimaterecovery3, Lshape, crack
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

[Dlambda, area] = gradbasis(node,elem);
gradU = gradu(node,elem,u,Dlambda); % Save the piecewise constant gradient
Du = recovery(node,elem,gradU,area); % Recover to piecewise linear

% Compute the relative error estimator eta_r if requested
if nargout >= 3
    nt = size(elem,1);
    numerator = 0; % ||R_h(u_h^n) - ∇u_h^n||^2
    denominator = 0; % ||R_h(u_h^n)||^2
    
    for i = 1:nt
        % Indices of the three nodes of the current element
        idx = elem(i,:);
        aT = area(i); % Area of the current element
        
        % Recovered gradient at the three nodes (3x2 matrix: [du/dx, du/dy])
        DuNodes = Du(idx,:);
        
        % Original piecewise constant gradient on the element (1x2 vector)
        gradU_elem = gradU(i,:);
        
        % Extract x and y components
        Du_x = DuNodes(:,1); % x-component at three nodes
        Du_y = DuNodes(:,2); % y-component at three nodes
        
        % Sum and sum of squares for x-component
        Sx = sum(Du_x);
        Qx = sum(Du_x.^2);
        % Integral of (R_h,x - gradU_x)^2 over the element
        int_x = aT * ( (Qx + Sx^2)/12 - (2/9)*gradU_elem(1)*Sx + gradU_elem(1)^2 );
        
        % Sum and sum of squares for y-component
        Sy = sum(Du_y);
        Qy = sum(Du_y.^2);
        % Integral of (R_h,y - gradU_y)^2 over the element
        int_y = aT * ( (Qy + Sy^2)/12 - (2/9)*gradU_elem(2)*Sy + gradU_elem(2)^2 );
        
        % Add to numerator
        numerator = numerator + int_x + int_y;
        
        % Integral of |R_h|^2 = (R_h,x)^2 + (R_h,y)^2 over the element
        int_denom_x = aT * (Qx + Sx^2) / 12;
        int_denom_y = aT * (Qy + Sy^2) / 12;
        denominator = denominator + int_denom_x + int_denom_y;
    end
    
    % Avoid division by zero
    if denominator < eps
        eta_r = 0;
    else
        eta_r = sqrt(numerator) / sqrt(denominator);
    end
end

% Original estimator based on second derivatives
DDu = zeros(size(elem,1),4);
DDu(:,1:2) = gradu(node,elem,Du(:,1),Dlambda);
DDu(:,3:4) = gradu(node,elem,Du(:,2),Dlambda);
eta = area.*sum(abs(DDu),2);
end