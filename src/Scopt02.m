% 求解器主入口函数
%
% 输入
% xi : 
% 
% 输出
% Xi        : 
% invXi     : 
% 

function [ x, y, s ] = Scopt02( A, b, c, K )

% 设置求解器参数
maxIter = 500 ;             % 最大迭代步数
epsilon = 1e-9  ;           % 对偶间隙误差限
gamma   = 0.995 ;           % 固定常数
bigNum  = 1e10  ;           % 固定大常数


if ~isfield( K, 'l' ), K.l = 0; end     % 判别线性锥是否存在
if ~isfield( K, 'q' ), K.q = 0; end     % 判别二阶锥是否存在
if ~isfield( K, 'r' ), K.r = 0; end     % 判别旋转二阶锥是否存在

% ==========================
% 初始化赋值
% ==========================
m  = size( b, 1 ) ;                     % 对偶决策变量 y 的维数
n  = K.l + sum( K.q ) + sum( K.r ) ;    % 原始决策变量 x 的维数

e1 = zeros( n, 1 ) ;                    % 单位列向量 e
T  = zeros( n, n ) ;                    % 转换矩阵 T
Q  = zeros( n, n ) ;                    % 转换矩阵 Q
x  = zeros( n, 1 ) ;                    % 原始决策变量 x in K
s  = zeros( n, 1 ) ;                    % 对偶松弛变量 s in K

if K.l > 0                              % 线性锥情况
    for k = 1: K.l
        i = k ;                         % 行列标
        [ e1( i ), T( i, i ), Q( i, i ) ] = TransMatLinearCone( k ) ;
        x( i ) = T( i, i )*e1( i ) ;
        s( i ) = T( i, i )*e1( i ) ;
    end
end

if K.q( 1 ) > 0                         % 二阶锥情况
    for k = 1: length( K.q )
        i = [ K.l + sum( K.q( 1: k ) ) - K.q( k ) + 1: ...
                K.l + sum( K.q( 1: k ) ) ] ;
        [ e1( i ), T( i, i ), Q( i, i ) ] = TransMatQuadCone( K.q( k ) ) ;
        x(  i ) = T( i, i )*e1( i ) ;
        s(  i ) = T( i, i )*e1( i ) ;
    end
end

if K.r( 1 ) > 0                         % 旋转二阶锥情况
    for k = 1: length( K.r )
        i = [ K.l + sum( K.q ) + sum( K.r( 1: k ) ) - K.r( k ) + 1: ...
                K.l + sum( K.q ) + sum( K.r( 1: k ) ) ] ;
        [ e1( i ), T( i, i ), Q( i, i ) ] = TransMatRotQuadCone( K.r( k ) ) ;
        x(  i ) = T( i, i )*e1( i ) ;
        s(  i ) = T( i, i )*e1( i ) ;
    end
end

Deg   = e1'*e1 ;         % 对称锥的度数 Deg
y     = zeros( m, 1 ) ;  % 对偶决策变量 y
tau   = 1 ;              % Goldman-Tucker齐次模型附加变量 tau
kappa = 1 ;              % 对偶间隙松弛变量 kappa

iter  = 1 ;              % 迭代步数计数器
alpha = 0 ;

tic ;
% while x'*s + tau*kappa > epsilon
for iter = 1: maxIter
    mu0 = ( x'*s + tau*kappa )/( Deg + 1 ) ; % 中心参数初始值
    
    if ( mod(iter, 2) ~= 0 )
        % ===========================
        % 预估步
        % ===========================
        nu        = 0 ;
        Exs       = zeros( n, 1 ) ;
        Ekappatau = 0 ;
        g         = x'*s + kappa*tau ;
    else
        % ===========================
        % 校正步
        % ===========================
        gp = ( x + alpha*dx )'*( s + alpha*ds ) ...
           + ( kappa + alpha*dkappa )*( tau + alpha*dtau ) ;
        nu = ( ( gp )/g )^2*( gp/( Deg + 1 ) ) ;
%         nu = ( gp/g )^3
        
        dXTilde   = zeros( n, n ) ;             % 决策向量 x 牛顿方向 dx 的矩阵化
        dSTilde   = zeros( n, n ) ;             % 决策向量 s 牛顿方向 ds 的矩阵化
        if K.l > 0                              % 线性锥情况
            for k = 1: K.l
                i = k ;
                dXTilde( i, i ) = x( i ) ;      % dXTilde = mat( T*( Theta*G )*dx )
                dSTilde( i, i ) = s( i ) ;      % dSTilde = mat( T*( Theta*G )^( -1 )*ds )
            end
        end

        if K.q(1) > 0                           % 二阶锥情况
            for k = 1: length( K.q )
                i = [ K.l + sum( K.q( 1: k ) ) - K.q( k ) + 1: ...
                      K.l + sum( K.q( 1: k ) ) ] ;
                [ dXTilde( i, i ), ~ ] = ArrowHeadMat( ...
                    T( i, i )*( Theta( i, i )*G( i, i ) )*dx( i ) ) ;
                [ dSTilde( i, i ), ~ ] = ArrowHeadMat( ...
                    T( i, i )*( invG( i, i )*invTheta( i, i ) )*ds( i ) ) ;
            end
        end

        if K.r(1) > 0                           % 旋转二阶锥情况
            for k = 1: length( K.r )
                i = [ K.l + sum( K.q ) + sum( K.r( 1: k ) ) - K.r( k ) + 1: ...
                      K.l + sum( K.q ) + sum( K.r( 1: k ) ) ] ;
                [ dXTilde( i, i ), ~ ] = ArrowHeadMat( ...
                    T( i, i )*( Theta( i, i )*G( i, i ) )*dx( i ) ) ;
                [ dSTilde( i, i ), ~ ] = ArrowHeadMat( ...
                    T( i, i )*( invG( i, i )*invTheta( i, i ) )*ds( i ) ) ;
            end
        end
    
        Exs       = dXTilde*dSTilde*e1 ;         % 修正项
        Ekappatau = dkappa*dtau ;                % 修正项
%         if abs(dd) < epsilon
%             dtau = 0 ;
%         else
%             dtau=(rrt-b2'*dy)/dd;
%         end;
    end
    
    % ===============================
    % Nesterov-Todd 放缩
    % ===============================
    Theta     = zeros( n, n ) ;             % 正放缩系数矩阵
    invTheta  = zeros( n, n ) ;             % 正放缩系数矩阵的逆矩阵
    G         = zeros( n, n ) ;             % 放缩矩阵
    invG      = zeros( n, n ) ;             % 放缩矩阵的逆矩阵
    XTilde    = zeros( n, n ) ;             % 决策向量 x 的矩阵化
    invXTilde = zeros( n, n ) ;             % 矩阵 XTilde 的逆矩阵
    STilde    = zeros( n, n ) ;             % 决策向量 s 的矩阵化
    if K.l > 0                              % 线性锥情况
        for k = 1: K.l
            i = k ;
            Theta(    i, i ) = 1 ;
            invTheta( i, i ) = 1 ;
            G(        i, i ) = 1 ;
            invG(     i, i ) = 1 ;
            XTilde(   i, i )    = x( i ) ;     % XTilde = mat( T*( Theta*G )*x )
            invXTilde(   i, i ) = 1/x( i ) ;
            STilde(   i, i )    = s( i ) ;     % STilde = mat( T*( Theta*G )^( -1 )*s )
        end
    end
    
    if K.q(1) > 0                           % 二阶锥情况
        for k = 1: length( K.q )
            i = [ K.l + sum( K.q( 1: k ) ) - K.q( k ) + 1: ...
                  K.l + sum( K.q( 1: k ) ) ] ;
            [ Theta( i, i ), invTheta( i, i ), ...
                G( i, i ), invG( i, i ) ] = ScaleMatQuadCone( ...
                e1( i ), x( i ), s( i ), ...
                T( i, i ), Q( i, i ), K.q( k ) ) ;
            [ XTilde( i, i ), invXTilde( i, i ) ] = ArrowHeadMat( ...
                T( i, i )*( Theta( i, i )*G( i, i ) )*x( i ) ) ;
            [ STilde( i, i ), ~ ] = ArrowHeadMat( ...
                T( i, i )*( invG( i, i )*invTheta( i, i ) )*s( i ) ) ;

        end
    end
    
    if K.r(1) > 0                           % 旋转二阶锥情况
        for k = 1: length( K.r )
            i = [ K.l + sum( K.q ) + sum( K.r( 1: k ) ) - K.r( k ) + 1: ...
                  K.l + sum( K.q ) + sum( K.r( 1: k ) ) ] ;
            [ Theta( i, i ), invTheta( i, i ), ...
                G( i, i ), invG( i, i ) ] = ScaleMatRotQuadCone( ...
                e1( i ), x( i ), s( i ), ...
                T( i, i ), Q( i, i ), K.r( k ) ) ;
            [ XTilde( i, i ), invXTilde( i, i ) ] = ArrowHeadMat( ...
                T( i, i )*( Theta( i, i )*G( i, i ) )*x( i ) ) ;
            [ STilde( i, i ), ~ ] = ArrowHeadMat( ...
                T( i, i )*( invG( i, i )*invTheta( i, i ) )*s( i ) ) ;

        end
    end
    
    rp =  A*x - b*tau ;              % 原始不可行性测量
    rd = -A'*y + c*tau - s ;        % 对偶不可行性测量
    rg =  b'*y - c'*x - kappa ;      % 互补间隙测量

    r1 = rp*nu - (  A*x - b*tau ) ;
    r2 = rd*nu - ( -A'*y + c*tau - s ) ;
    r3 = rg*nu - (  b'*y - c'*x - kappa ) ;
    r4 = nu*mu0*e1 - XTilde*STilde*e1 - Exs ;
    r5 = nu*mu0 - kappa*tau - Ekappatau ;
    
    % ==================================
    % LU分解
    % ==================================
    M = [ A                   , -b           , zeros( m, m )  , zeros( m, n )             , zeros( m, 1 ) ;
          zeros( n, n )       , c            , -A'            , -eye( n, n )              , zeros( n, 1 ) ;
          -c'                 , 0            , b'             , zeros( 1, n )             , -1            ;
          STilde*T*( Theta*G ), zeros( n, 1 ), zeros( n, m )  , XTilde*T*( invG*invTheta ), zeros( n, 1 ) ;
          zeros( 1, n )       , kappa        , zeros( 1, m )  , zeros( 1, n )             , tau           ; ] ;
    F = [ r1 ; r2 ; r3 ; r4 ; r5 ] ;

    [ L, U ] = lu( M ) ;        % LU 分解求线性方程组
    Delta = U\( L\F ) ;
    
    dx     = Delta( 1: n ) ;
    dtau   = Delta( n + 1 ) ;
    dy     = Delta( n + 2: n + m + 1 ) ;
    ds     = Delta( n + m + 2: 2*n + m + 1 ) ;
    dkappa = Delta( 2*n + m + 2 ) ;

    % ================================
    % 计算牛顿步长 alpha
    % ================================
    % dx
    alphax = 1.0 ;
    if K.l > 0                              % 线性锥情况
        for k = 1: K.l
            i = k ;
            if dx( i ) >= 0
                alphax = min( [ alphax, bigNum ] ) ;
            else
                alphax = min( [ alphax, -x( i )/dx( i ) ] ) ;
            end
        end
    end

    if K.q( 1 ) > 0                           % 二阶锥情况
        for k = 1: length( K.q )
            i = [ K.l + sum( K.q( 1: k ) ) - K.q( k ) + 1: ...
                  K.l + sum( K.q( 1: k ) ) ] ;
            gamma1 = dx( i( 1 ) )^2 ...
                   - dx( i( 2: K.q( k ) ) )'*dx( i( 2: K.q( k ) ) ) ;
            gamma2 = 2*( x( i( 1 ) )*dx( i( 1 ) ) ...
                   - x( i( 2: K.q( k ) ) )'*dx( i( 2: K.q( k ) ) ) ) ;
            gamma3 = x( i( 1 ) )^2 ...
                   - x( i( 2: K.q( k ) ) )'*x( i( 2: K.q( k ) ) ) ;
            if gamma2^2 - 4*gamma1*gamma3 < 0
                lambdaTilde = bigNum ;
            else
                if gamma1 > 0
                    lambda1 = ( -gamma2 + sqrt( gamma2^2 - 4*gamma1*gamma3 ) )/( 2*gamma1 ) ;
                    lambda2 = ( -gamma2 - sqrt( gamma2^2 - 4*gamma1*gamma3 ) )/( 2*gamma1 ) ;
                    lambdaTilde = min( [ lambda1, lambda2 ] ) ;
                elseif gamma1 == 0
                    if gamma2 >= 0
                        lambdaTilde = bigNum ;
                    else
                        lambdaTilde = -gamma3/gamma2 ;
                    end
                else
                    lambda1 = ( -gamma2 + sqrt( gamma2^2 - 4*gamma1*gamma3 ) )/( 2*gamma1 ) ;
                    lambda2 = ( -gamma2 - sqrt( gamma2^2 - 4*gamma1*gamma3 ) )/( 2*gamma1 ) ;
                    lambdaTilde = max( [ lambda1, lambda2 ] ) ;
                end
            end
            if dx( i( 1 ) ) >= 0
                lambdaBar = bigNum ;
            else
                lambdaBar = -x( i( 1 ) )/dx( i( 1 ) ) ;
            end
            alphax = min( [ alphax, lambdaTilde, lambdaBar ] ) ;
        end
    end
    % ds
    alphas = 1.0 ;
    if K.l > 0                              % 线性锥情况
        for k = 1: K.l
            i = k ;
            if ds( i ) >= 0
                alphas = min( [ alphas, bigNum ] ) ;
            else
                alphas = min( [ alphas, -s( i )/ds( i ) ] ) ;
            end
        end
    end

    if K.q( 1 ) > 0                           % 二阶锥情况
        for k = 1: length( K.q )
            i = [ K.l + sum( K.q( 1: k ) ) - K.q( k ) + 1: ...
                  K.l + sum( K.q( 1: k ) ) ] ;
            gamma1 = ds( i( 1 ) )^2 ...
                   - ds( i( 2: K.q( k ) ) )'*ds( i( 2: K.q( k ) ) ) ;
            gamma2 = 2*( s( i( 1 ) )*ds( i( 1 ) ) ...
                   - s( i( 2: K.q( k ) ) )'*ds( i( 2: K.q( k ) ) ) ) ;
            gamma3 = s( i( 1 ) )^2 ...
                   - s( i( 2: K.q( k ) ) )'*s( i( 2: K.q( k ) ) ) ;
            if gamma2^2 - 4*gamma1*gamma3 < 0
                lambdaTilde = bigNum ;
            else
                if gamma1 > 0
                    lambda1 = ( -gamma2 + sqrt( gamma2^2 - 4*gamma1*gamma3 ) )/( 2*gamma1 ) ;
                    lambda2 = ( -gamma2 - sqrt( gamma2^2 - 4*gamma1*gamma3 ) )/( 2*gamma1 ) ;
                    lambdaTilde = min( [ lambda1, lambda2 ] ) ;
                elseif gamma1 == 0
                    if gamma2 >= 0
                        lambdaTilde = bigNum ;
                    else
                        lambdaTilde = -gamma3/gamma2 ;
                    end
                else
                    lambda1 = ( -gamma2 + sqrt( gamma2^2 - 4*gamma1*gamma3 ) )/( 2*gamma1 ) ;
                    lambda2 = ( -gamma2 - sqrt( gamma2^2 - 4*gamma1*gamma3 ) )/( 2*gamma1 ) ;
                    lambdaTilde = max( [ lambda1, lambda2 ] ) ;
                end
            end
            if ds( i( 1 ) ) >= 0
                lambdaBar = bigNum ;
            else
                lambdaBar = -s( i( 1 ) )/ds( i( 1 ) ) ;
            end
            alphas = min( [ alphas, lambdaTilde, lambdaBar ] ) ;
        end
    end
%     if K.r(1) > 0                           % 旋转二阶锥情况
%         for k = 1: length( K.r )
%             i = [ K.l + sum( K.q ) + sum( K.r( 1: k ) ) - K.r( k ) + 1: ...
%                   K.l + sum( K.q ) + sum( K.r( 1: k ) ) ] ;
%             
%         end
%     end
    
    if dtau >= 0
        alphatau = bigNum ;
    else
        alphatau = -tau/dtau ;
    end
    
    if dkappa >= 0
        alphakappa = bigNum ;
    else
        alphakappa = -kappa/dkappa ;
    end
    
    alpha = min( [ 1.0              ; ...
                   gamma*alphax     ; ...
                   gamma*alphas     ; ...
                   gamma*alphatau   ; ...
                   gamma*alphakappa ; ...
                   ] ) ;
    
    % ================================
    % 更新步
    % ================================
    if ( mod(iter, 2) <= 0 )        % 校正步更新
        x     = x     + alpha*dx     ;
        tau   = tau   + alpha*dtau   ;
        y     = y     + alpha*dy     ;
        s     = s     + alpha*ds     ;
        kappa = kappa + alpha*dkappa ;
    end
%     iter = iter + 1 ;
    
    % 最大迭代步数
    if iter >= maxIter
        x = x/tau ;
        y = y/tau ;
        s = x/tau ;
        fprintf('maximum iteration is reached ! iter = %d\n', floor( iter/2 ) ) ;
        break ;
    end
    if x'*s + tau*kappa <= epsilon
        break ;
    end
end
toc

% ================================
% 解分析
% ================================
if  tau > 0
    x = x/tau ;
    y = y/tau ;
    s = x/tau ;
    fprintf( 'primal-dual optimal solution is found !\n' ) ;
    fprintf( 'iteration number is %d\n', floor( iter/2 ) ) ;
elseif kappa == 0
    fprintf( 'problem is ill-posed !\n' ) ;
elseif c'*x < 0
    fprintf( 'dual problem is infeasible !\n' ) ;
elseif b'*y > 0
    fprintf( 'primal problem is infeasible !\n' ) ;
else
    fprintf( 'The problem is either infeasible or unbounded !\n' ) ;
end



end



