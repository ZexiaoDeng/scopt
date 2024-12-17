%  设置第 i 个线性锥的 NT 放缩系数 theta 和放缩矩阵 G
% 输入
% k : 锥维数
% 输出
% theta : 放缩系数
% G     : 放缩矩阵
% 
function [ Thetai, invThetai, Gi, invGi ] = ScaleMatLinearCone( k )

if k > 0
    thetai    = 1 ;
    Thetai    = 1 ;
    invThetai = 1 ;
    Gi        = 1 ;
    invGi     = 1 ;
else
    Thetai    = [] ;
    invThetai = [] ;
    Gi        = [] ;
    invGi     = [] ;
end
