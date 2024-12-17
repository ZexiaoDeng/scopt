
clc
clear
close all

xi = [ sqrt( 2 ), 1, 0 ]' ;
si = [ sqrt( 2 ), 0, 1 ]' ;

for k = 1 : 3
    [ Thetai, invThetai, Gi, invGi ] = ScaleMatLinearCone( k )
end
Gi*invGi

k = 3 ;
[ e1i, Ti, Qi ] = TransMatQuadCone( k )
[ Thetai, invThetai, Gi, invGi ] = ScaleMatQuadCone( e1i, xi, si, Ti, Qi, k )
Thetai*invThetai
Gi*invGi

% 上面的案例 si 并不在旋转锥二阶锥内，需要修改一下
xi = [ sqrt( 2 ), 1, 0 ]' ;
si = [ sqrt( 2 ), 1, 1 ]' ;
k = 3 ;
[ e1i, Ti, Qi ] = TransMatRotQuadCone( k )
[ Thetai, invThetai, Gi, invGi ] = ScaleMatRotQuadCone( e1i, xi, si, Ti, Qi, k )
Thetai*invThetai
Gi*invGi

