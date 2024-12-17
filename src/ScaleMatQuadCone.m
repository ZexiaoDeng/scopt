%  设置第 i 个二阶锥的 NT 放缩系数 theta 和放缩矩阵 G
% 输入
% k : 锥维数
% 输出
% theta : 放缩系数
% G     : 放缩矩阵
% 
function [ Thetai, invThetai, Gi, invGi ] = ScaleMatQuadCone( e1i, xi, si, Ti, Qi, k )

if k > 0 & k < 2
    error( 'Input error for quadratic cone!' ) ;
elseif k >= 2
    sthetai   = sqrt( ( si'*Qi*si )/( xi'*Qi*xi ) ) ;
    thetai    = sqrt( sthetai ) ;
    Thetai    = thetai*eye( k ) ;
    invThetai = ( 1/thetai )*eye( k ) ;
    
    gi        = ( ( 1/thetai )*si + thetai*Qi*xi ) ...
              / ( sqrt( 2 )*sqrt( xi'*si + sqrt( xi'*Qi*xi*si'*Qi*si ) ) ) ;
    Gi        = -Qi + ( e1i + gi )*( e1i + gi )'/( 1 + e1i'*gi ) ;
    invGi     = Qi*Gi*Qi ;
    
%     % 验证
%     invGi  = Qi*Gi*Qi ;
%     sGi    = -Qi + 2*gi*gi' ;
%     sinvGi = -Qi + 2*( Qi*gi )*( Qi*gi )' ;
%     si     = thetai^2*Gi^2*xi 
%     
%     Gi*invGi
%     norm( Gi*Gi - sGi )
%     norm( Gi*Gi*invsGi )
else
    Thetai    = [] ;
    invThetai = [] ;
    Gi        = [] ;
    invGi     = [] ;
end