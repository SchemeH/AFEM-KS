function KS_2D_afem()


% 参数设置
global alpha;
alpha = 1;

% 域和分区
xb = 0; xe = 1;
yb = 0; ye = 1;
T0 = 0; Te = 1e-4; 
N = 32; M = N; NK = 10;
hx = (xe-xb)/N; hy = (ye-yb)/M; dt = (Te-T0)/NK; TK = NK;
T = T0:dt:Te;
% 计算资源限制
maxdof = 512*512; 

% 自适应参数
Tol_Space = 1e-1; 
theta = 0.7;       % 细化参数
ctheta = 0.4;      % 粗化参数   

fprintf(1,'\n *************************************************\n');
fprintf(1,'\n --- Parameters and Mesh Sizes ---\n');
fprintf(1,'\n Tolerence of space = %e, Maximum Dof = %e, dt = %e\n',Tol_Space, maxdof, dt);
fprintf(1,'\n theta = %e, ctheta = %e, alpha = %e\n',theta,ctheta,alpha);
fprintf(1,'\n hx = %d, hy = %d, dt = %d\n',hx,hy,dt);


% 初始网格和边界
[node,elem] = squaremesh([xb xe yb ye], min(hx,hy));
bdFlag = setboundary(node,elem,'Neumann');
pde = KSdata2d;

% 初始网格细化
n = 0;
while true
    uold = pde.u0(node);
    [eta, ~, ~] = relativeerror3(node,elem,uold,'max');
    TotalError = norm(eta,inf);
    Dof = size(node, 1);
    fprintf(1,'\n refine times = %d, dof = %d, error_space = %e\n',n,Dof,TotalError);
    if TotalError < 1e-1 || Dof > maxdof
        break;
    else
        markedElem = mark(elem,eta,theta,'MAX');
        n = n + 1;
    end
    [node,elem,bdFlag,HB] = bisect(node,elem,markedElem,bdFlag);
end

fprintf(1,'\n *************************************************\n');
% 初始值
uold = pde.u0(node); vold = pde.v0(node);

% visualization
% refine mesh
figure(1);
showmesh(node,elem)
title('Initial mesh')
title('t=0')

% initial values on refine mesh
figure(2);
trisurf(elem,node(:,1),node(:,2),uold);
shading interp
colormap jet
colorbar;
xlabel('X');
ylabel('Y');
title('Initial U')

figure(3);
trisurf(elem,node(:,1),node(:,2),vold);
shading interp
colormap jet
colorbar;
xlabel('X');
ylabel('Y');
title('Initial V');

% 预分配存储
dofs = 0:1:TK; 
uinf = 0:1:TK; vinf = 0:1:TK; 
umin = 0:1:TK; vmin = 0:1:TK; 
total_mass = 0:1:TK; egy = 0:1:TK; 

dofs(1) = size(node, 1);
uinf(1) = max(uold); vinf(1) = max(vold); 
umin(1) = min(uold); vmin(1) = min(vold);
total_mass(1) = compute_element_mass_fem(node, elem, uold);
egy(1) = get_energy(node, elem, uold, vold, alpha);

t0 = cputime;
for i = 1:NK
    % 自适应时间步长
    [u1, v1, fu, fv, M, A, Ac] = KS_onestep_ETD_first(elem, node, i*dt, dt, uold, vold, pde, pde);
    u1 = nonnegative_projection(u1, total_mass(1), node, elem, 1e-8); v1 = max(0,v1);
    [u, v] = KS_onestep_ETD_second(elem, node, dt, uold, vold, u1, v1, fu, fv, M, A, Ac);
    u = nonnegative_projection(u, total_mass(1), node, elem, 1e-8);  v = max(0,v);

    % 空间自适应模块 - 使用增强版
    [node, elem, u, v, HB, bdFlag] = spaceAdaptation_enhanced(...
        node, elem, u, v, bdFlag, HB,...
        'TolSpace', 0.1, ...
        'MaxDof', maxdof, ...
        'RefineTheta', theta, ...
        'CoarsenTheta', ctheta);

    % 更新解
    u = nonnegative_projection(u, total_mass(1), node, elem, 1e-8); v = max(0,v);
    uold = u; vold = v;
    
    % 记录变量
    dofs(i + 1) = size(node, 1); uinf(i + 1) = max(u);  vinf(i + 1) = max(v);
    umin(i + 1) = min(u);  vmin(i + 1) = min(v);
    total_mass(i+1)= compute_element_mass_fem(node, elem, uold);
    egy(i+1) = get_energy(node, elem, u, v, alpha);     
end
t1 = cputime - t0;
% Output the recorded values
fprintf(1,'\n *************************************************\n');
fprintf (1,'\n MY_PROGRAM took %f seconds to run !!!\n',t1);
fprintf(1,'\n Minimal value of U = %f, Maximal value of U = %f\n',max(uinf),min(umin));
fprintf(1,'\n Minimal value of V = %f, Maximal value of V = %f\n',max(vinf),min(vmin));
fprintf(1,'\n Initial energy = %e, Final energy is %e\n',egy(1),egy(end));
fprintf(1,'\n Initial volume of U = %e, Final volume of U is %e\n',total_mass(1),total_mass(end));


% visualization
figure(4);
showmesh(node,elem);
t1 = ['t=' , num2str(i*dt) , ',dof=' , num2str(size(node, 1))];
title(t1);

figure(5);
trisurf(elem,node(:,1),node(:,2),u);
shading interp
colormap jet
colorbar;
xlabel('X');
ylabel('Y');
t2 = ['t=' , num2str(i*dt), ',Uh'];
title(t2);

figure(6);
trisurf(elem,node(:,1),node(:,2),v);
shading interp
colormap jet
colorbar;
xlabel('X');
ylabel('Y');
t2 = ['t=' , num2str(i*dt), ',Vh'];
title(t2);

figure(7);
plot(T,uinf,'.-');
ylabel('Maximum norm: U')
xlabel('Time');

figure(8);
plot(T,vinf,'.-');
ylabel('Maximum norm: V')
xlabel('Time');

figure(9);
plot(T,vmin,'.-');
hold on;
plot(T,umin,'-');
ylabel('Minimum norm')
xlabel('Time');
legend('minimum values of V','minimum values of U');

figure(10);
plot(T,total_mass-total_mass(1),'.-');
ylabel('Mass error')
xlabel('Time');

figure(11);
plot(T,egy,'.-');
ylabel('Energy')
xlabel('Time');


figure(12);
plot(T,dofs,'.-');
ylabel('adaptive spatial mesh')
xlabel('Time');
end

% ============= 增强版空间自适应函数 =============
function [node, elem, u, v, HB, bdFlag] = spaceAdaptation_enhanced(...
    node, elem, u, v, bdFlag, HB, varargin)
% 参数解析
p = inputParser;
addParameter(p, 'TolSpace', 1e-3);     % 空间误差容限
addParameter(p, 'MaxDof', 1e4);        % 最大自由度限制
addParameter(p, 'RefineTheta', 0.7);   % 细化阈值
addParameter(p, 'CoarsenTheta', 0.3);  % 粗化阈值
parse(p, varargin{:});

% 初始化变量
iter_count = 0;
prev_dof = size(node, 1);
converged = false;

% 主自适应循环 - 移除迭代次数限制
while ~converged
    iter_count = iter_count + 1;
    
    % 1. 计算当前空间误差
    % [eta, ~, TotalError] = relativeerror3(node,elem,u,'max');
    [eta, ~, ~] = relativeerror3(node,elem,uold,'max');
    TotalError = norm(eta,inf);
    Dof = size(node, 1);
    
    fprintf('\nAdapt Iter %d: Dof=%d, Error=%.4e (Target=%.4e)', ...
            iter_count, Dof, TotalError, p.Results.TolSpace);
    
    % 2. 决策逻辑 - 基于容差的细化/粗化
    if TotalError > p.Results.TolSpace
        % 误差大于容差 -> 执行细化
        markedElem = mark(elem, eta, p.Results.RefineTheta, 'MAX');
        
        if ~isempty(markedElem)
            [node, elem, bdFlag, HB] = bisect(node, elem, markedElem, bdFlag);
            u = nodeinterpolate(u, HB);
            v = nodeinterpolate(v, HB);
            fprintf(' - Refined %d elements', numel(markedElem));
        else
            fprintf(' - No elements to refine');
            converged = true; % 无法继续细化
        end
        
    elseif TotalError < p.Results.TolSpace
        % 误差显著小于容差 -> 执行粗化
        [markedElem, canCoarsen] = mark_coarsen(elem, eta, p.Results.CoarsenTheta);
        
        if canCoarsen
            [node, elem, bdFlag, indexMap] = coarsen(node, elem, markedElem, bdFlag);
            u = nodeinterpolate(u, indexMap);
            v = nodeinterpolate(v, indexMap);
            fprintf(' - Coarsened %d elements', numel(markedElem));
        else
            fprintf(' - Cannot coarsen further');
            converged = true; % 无法继续粗化
        end
        
    else
        % 误差在容差附近 (±20%) -> 保持当前网格
        fprintf(' - Error within tolerance range');
        converged = true;
    end
    
    % 3. 检查收敛条件
    current_dof = size(node, 1);
    if converged
        fprintf(' - Adaptation converged');
        break;
    end
    
    % 4. 检查最大自由度限制
    if current_dof >= p.Results.MaxDof
        fprintf(' - Reached maximum DOF limit (%d)', p.Results.MaxDof);
        break;
    end
    
    % 5. 检查网格变化
    if current_dof == prev_dof
        fprintf(' - Mesh unchanged');
        converged = true;
    else
        prev_dof = current_dof;
    end
end

fprintf('\nFinal: Dof=%d, Error=%.4e, Iterations=%d\n', ...
        size(node, 1), TotalError, iter_count);
end

function [markedElem, canCoarsen] = mark_coarsen(elem, eta, coarsenTheta)
% 标记可粗化元素，并返回是否能粗化
NT = size(elem, 1);
isMark = false(NT, 1);

% 策略1：基于阈值的标记
threshold = coarsenTheta * max(eta);
isMark(eta < threshold) = true;

% 策略2：确保至少标记20%的元素
minMarkCount = max(1, round(0.2 * NT));
if nnz(isMark) < minMarkCount
    [~, idx] = sort(eta, 'ascend');
    isMark(idx(1:minMarkCount)) = true;
end

% 策略3：检查是否能安全粗化（避免网格退化）
canCoarsen = true;
if nnz(isMark) < 2
    canCoarsen = false;
end

markedElem = uint32(find(isMark));
fprintf(' [Coarsen] Marked %.1f%% elements', 100*nnz(isMark)/NT);
end

function mass = compute_element_mass_fem(nodes, elements, u_values)
[~, M, ~] = assemblematrix(nodes, elements, 1); 
mass = u_values'*M;  
end

function egy = get_energy(nodes, elements, u, v, alpha)
[A, M, ~] = assemblematrix(nodes, elements, 1);  
u_safe = max(u, eps);  
f = u .* log(u_safe) - u - u .* v;          
integral_f = f'*M;                
integral_grad = (1/2) * alpha * (v' * A * v); 
egy = integral_f + integral_grad;           
end

function u_proj = nonnegative_projection(v, mass_target, node, elem, tol)
    % 初始化Lagrange乘子
    xi_prev = 0;
    xi_curr = 0;
    max_iter = 100;
    residual_prev = 0;
    
    for iter = 1:max_iter
        % 应用Lagrange乘子调整
        val = v - xi_curr;
        
        % 仅应用非负性约束（无上界约束）
        u_proj = val;
        u_proj(val < 0) = 0;  % 非负性约束
        
        % 计算当前总质量（使用质量加权点积）
        current_mass = compute_element_mass_fem(node, elem, u_proj);
        residual = current_mass - mass_target;
        
        % 收敛判断
        if abs(residual) < tol
            fprintf('投影收敛于迭代 %d, 残差 = %.6e\n', iter, residual);
            break;
        end
        
        % 割线法更新xi (Secant method)
        if iter == 1
            % 初始步长
            xi_next = xi_curr - residual * 0.1;
        else
            delta_xi = xi_curr - xi_prev;
            delta_res = residual - residual_prev;
            
            if abs(delta_res) < 1e-12
                fprintf('残差变化太小，终止于迭代 %d\n', iter);
                break;
            end
            
            % 割线法更新
            xi_next = xi_curr - residual * (delta_xi / delta_res);
        end
        
        % 更新迭代变量
        xi_prev = xi_curr;
        xi_curr = xi_next;
        residual_prev = residual;
        
        % 每10次迭代显示进度
        if mod(iter, 10) == 0
            fprintf('迭代 %d: xi = %.6e, 残差 = %.6e\n', iter, xi_curr, residual);
        end
    end
    
    if iter == max_iter
        warning('达到最大迭代次数 (%d)，未收敛', max_iter);
    end
    
    % 最终质量验证
    final_mass = compute_element_mass_fem(node, elem, u_proj);
    non_negative_violation = sum(u_proj < 0);
    
    fprintf('投影结果: 质量 = %.10e (目标 = %.10e), 差值 = %.10e, 负值节点数 = %d\n',...
            final_mass, mass_target, final_mass - mass_target, non_negative_violation);
end

