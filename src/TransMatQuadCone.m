%  设置第 i 个二阶锥的转换矩阵 T 和 Q
% 输入
% k : 锥维数
% 输出
% T : 转换矩阵
% Q : 转换矩阵
function [ e1i, Ti, Qi ] = TransMatQuadCone( k )

if k > 0 & k < 2
    error( 'Input error for quadratic cone!' ) ;
elseif k >= 2
    e1i = [ 1 ; zeros( k - 1, 1 ) ; ] ;
    Ti  =  eye( k ) ;
    Qi  = blkdiag( 1, -eye( k - 1 ) ) ;
else
    e1i = [] ;
    Ti  = [] ;
    Qi  = [] ;
end

end

