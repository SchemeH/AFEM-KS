function KS_2D_afem_v1()


% 参数设置
global alpha;
alpha = 1;

% 域和分区
xb = 0; xe = 1;
yb = 0; ye = 1;
T0 = 0; Te = 1e-3; 
N = 32; M = N; NK = 100;
hx = (xe-xb)/N; hy = (ye-yb)/M; dt = (Te-T0)/NK; 

% 计算资源限制
maxnk = 1e6; maxdof = 256*256; mindt = 1/maxnk; maxdt = 1e-3; 

% 自适应参数
Tol_Space = 0.1; 
theta = 0.7;       % 细化参数
ctheta = 0.4;      % 粗化参数   

fprintf(1,'\n *************************************************\n');
fprintf(1,'\n --- Parameters and Mesh Sizes ---\n');
fprintf(1,'\n Tolerence of space = %e\n',Tol_Space);
fprintf(1,'\n theta = %e, ctheta = %e, alpha = %e\n',theta,ctheta,alpha);
fprintf(1,'\n Maximum Dof = %e, Minimum time step size = %e, Maximum time step size = %e\n',maxdof,mindt,maxdt);
fprintf(1,'\n hx = %d, hy = %d, dt = %d\n',hx,hy,dt);
fprintf(1,'\n N = %d, M = %d, NK = %d\n',N,M,NK);

% 初始网格和边界
[node,elem] = squaremesh([xb xe yb ye], min(hx,hy));
bdFlag = setboundary(node,elem,'Neumann');
pde = KSdata2d;

% 初始网格细化
n = 0;
while true
    uold = pde.u0(node);
    [eta, ~, TotalError] = relativeerror3(node,elem,uold,'max');
    Dof = size(node, 1);
    fprintf(1,'\n refine times = %e, dof = %e, error_space = %e\n',n,Dof,TotalError);
    if TotalError < Tol_Space || Dof > maxdof
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

% % visualization
% % refine mesh
% figure(1);
% showmesh(node,elem)
% title('t=0')
% 
% % initial values on refine mesh
% figure(2);
% trisurf(elem,node(:,1),node(:,2),uold);
% shading interp
% colormap jet
% colorbar;
% xlabel('X');
% ylabel('Y');
% title('\rho_h');
% 
% figure(3);
% trisurf(elem,node(:,1),node(:,2),vold);
% shading interp
% colormap jet
% colorbar;
% xlabel('X');
% ylabel('Y');
% title('c_h');

% 预分配存储
dofs = zeros(1,maxnk); 
uinf = zeros(1,maxnk);  vinf = zeros(1,maxnk); 
umin = zeros(1,maxnk);  vmin = zeros(1,maxnk); 
relative_err = zeros(1,maxnk); adap_dt = zeros(1,maxnk);
Time_Error = zeros(1,maxnk); total_mass = zeros(1,maxnk);
egy = zeros(1,maxnk);

dofs(1) = size(node, 1);
uinf(1) = max(uold); vinf(1) = max(vold); 
umin(1) = min(uold); vmin(1) = min(vold);
relative_err(1) = 0;
adap_dt(1) = dt; Time_Error(1) = 0;
total_mass(1) = compute_element_mass_fem(node, elem, uold);
egy(1) = get_energy(node, elem, uold, vold, alpha);
egy_adap1(1) = egy(1); adap_dt(1) = dt;
tn_1 = 0; tn = 0;
i = 1; count = 1;
tic;
while tn<Te
    tn = tn_1 + dt;

    % 自适应时间步长
    [u1, v1, fu, fv, M, A, Ac] = KS_onestep_ETD_first(elem, node, tn, dt, uold, vold, pde, pde);
    u1 = nonnegative_projection(u1, total_mass(1), node, elem, eps);
    
    egy_adap1(i+1) = get_energy(node, elem, u1, v1, alpha);
    % 能量变化量与相对变化量
    delta_egy = egy_adap1(i+1) - egy_adap1(i);
    egy_grad = abs(delta_egy/adap_dt(i));
    De = abs(egy_grad)^2; 
    De = 10^6*De;
    dt = Adap_time_stepping(mindt, maxdt, De);

    [u, v] = KS_onestep_ETD_second(elem, node, dt, uold, vold, u1, v1, fu, fv, M, A, Ac);
    u = nonnegative_projection(u, total_mass(1), node, elem, eps);
    TimeError = max(max(abs(u-u1)))/max(max(abs(u)));
    
    % 空间自适应模块 - 使用增强版
    [node, elem, u, v, HB, bdFlag] = spaceAdaptation_enhanced(...
        node, elem, u, v, bdFlag, HB, total_mass(1), ...
        'TolSpace', Tol_Space, ...
        'MaxDof', maxdof, ...
        'RefineTheta', theta, ...
        'CoarsenTheta', ctheta);

    % 更新解
    uold = u; vold = v;
    u = nonnegative_projection(u, total_mass(1), node, elem, eps);
    % 记录变量
    dofs(i + 1) = size(node, 1); 
    uinf(i + 1) = max(u);  vinf(i + 1) = max(v);
    umin(i + 1) = min(u);  vmin(i + 1) = min(v);
    total_mass(i+1)= compute_element_mass_fem(node, elem, uold);
    egy(i+1) = get_energy(node, elem, u, v, alpha); 
    relative_err(i + 1) = TotalError;
    adap_dt(i + 1) = dt; 
    Time_Error(i + 1) = TimeError;
    
    i = i + 1; 
    tn_1 = tn; 
    count = count + 1;
    
    fprintf(1,'\n dtn = %e, dof = %e, Time = %e\n',dt,size(node, 1),tn);
end
t1 = toc;
% Post-processing
dofs = dofs(1:count); 
uinf = uinf(1:count);  vinf = vinf(1:count); 
umin = umin(1:count);  vmin = vmin(1:count); 
relative_err = relative_err(1:count); adap_dt = adap_dt(1:count);
Time_Error = Time_Error(1:count); egy = egy(1:count); 
total_mass = total_mass(1:count);
tstop = 0:(Te-T0)/(count-1):Te;

% Output the recorded values
fprintf(1,'\n *************************************************\n');
fprintf (1,'\n MY_PROGRAM took %f seconds to run !!!\n',t1);
fprintf(1,'\n Minimal value of U = %f, Maximal value of U = %f\n',max(uinf),min(umin));
fprintf(1,'\n Minimal value of V = %f, Maximal value of V = %f\n',max(vinf),min(vmin));
fprintf(1,'\n Initial energy = %e, Final energy is %e\n',egy(1),egy(count));
fprintf(1,'\n Initial volume of U = %e, Final volume of U is %e\n',total_mass(1),total_mass(count));


% visualization
figure(4);
tn = 1e-3;
showmesh(node,elem);
t1 = ['t=' , num2str(tn) , ',dof=' , num2str(size(node, 1))];
title(t1);

figure(5);
trisurf(elem,node(:,1),node(:,2),u);
shading interp
colormap jet
colorbar;
xlabel('X');
ylabel('Y');
title('\rho_h');

figure(6);
trisurf(elem,node(:,1),node(:,2),v);
shading interp
colormap jet
colorbar;
xlabel('X');
ylabel('Y');
title('c_h');

figure(7);
plot(tstop,uinf,'.-');
ylabel('Maximum norm: \rho_h')
xlabel('Time');

figure(8);
plot(tstop,vinf,'.-');
ylabel('Supremum norm: c_h')
xlabel('Time');

figure(9);
plot(tstop,vmin,'.-');
hold on;
plot(tstop,umin,'-');
ylabel('Minimum norm')
xlabel('Time');
legend('minimum values of \rho_h','minimum values of c_h');

figure(10);
plot(tstop,round(total_mass),'.-');
ylabel('Mass')
xlabel('Time');

figure(11);
plot(tstop,egy,'.-');
ylabel('Energy')
xlabel('Time');

figure(12);
plot(tstop,adap_dt,'.-');
ylabel('adaptive time step size')
xlabel('Time');

figure(13);
plot(tstop,dofs,'.-');
ylabel('adaptive spatial mesh')
xlabel('Time');

figure(14);
plot(tstop,Time_Error,'.-');
hold on;
plot(tstop,relative_err);
ylabel('Errors')
xlabel('Time');
legend('Time error','Space error');
end

% ============= 增强版空间自适应函数 =============
function [node, elem, u, v, HB, bdFlag] = spaceAdaptation_enhanced(...
    node, elem, u, v, bdFlag, HB, mass_target, varargin)
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
    [eta, ~, TotalError] = relativeerror3(node,elem,u,'max');
    Dof = size(node, 1);
    
    fprintf('\nAdapt Iter %d: Dof=%d, Error=%.4e (Target=%.4e)', ...
            iter_count, Dof, TotalError, p.Results.TolSpace);
    
    % 2. 决策逻辑 - 基于容差的细化/粗化
    if TotalError > p.Results.TolSpace
        % 误差大于容差 -> 执行细化
        markedElem = mark(elem, eta, p.Results.RefineTheta, 'MAX');
        
        if ~isempty(markedElem)
            [node, elem, bdFlag, HB] = bisect(node, elem, markedElem, bdFlag);
            u = nodeinterpolate(u, HB); u = nonnegative_projection(u, mass_target, node, elem, eps);
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
            u = nodeinterpolate(u, indexMap); u = nonnegative_projection(u, mass_target, node, elem, eps);
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
function V = Adap_time_stepping(dtmin,dtmax,De)
   v = dtmax/sqrt(1 + De);
   V = max(dtmin,v);
end

function u_proj = nonnegative_projection(v, mass_target, node, elem, tol)
    % 初始化搜索区间
    xi_low = min(v) - 1;    % 全正解（质量过大）
    xi_high = max(v);       % 全零解（质量=0）
    max_iter = 100;
    
    for iter = 1:max_iter
        xi = (xi_low + xi_high)/2;
        u_proj = max(v - xi, 0); % 平移+截断
        
        current_mass = compute_element_mass_fem(node, elem, u_proj);
        mass_diff = current_mass - mass_target;
        
        if abs(mass_diff) < tol
            break;
        end
        
        % 二分法更新区间
        if mass_diff > 0
            xi_low = xi;  % 质量过大→增大xi
        else
            xi_high = xi; % 质量过小→减小xi
        end
    end
    
    % 验证结果
    final_mass = compute_element_mass_fem(node, elem, u_proj);
    fprintf('投影: 目标质量=%.6e, 实际质量=%.6e, 差值=%.2e\n',...
            mass_target, final_mass, final_mass - mass_target);
end