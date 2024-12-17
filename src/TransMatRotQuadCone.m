%  设置第 i 个旋转二阶锥的转换矩阵 T 和 Q
% 输入
% k : 锥维数
% 输出
% T : 转换矩阵
% Q : 转换矩阵
function [ e1i, Ti, Qi ] = TransMatRotQuadCone( k )

if k > 0 & k < 3
    error( 'Input error for rotated quadratic cone!' ) ;
elseif k >= 3
    e1i = [ 1 ; zeros( k - 1, 1 ) ; ] ;
    Ti  =  blkdiag( [ 1/2^( 1/2 ),  1/2^( 1/2 ) ; ...
                     1/2^( 1/2 ), -1/2^( 1/2 ) ; ], ...
                     eye( k - 2 ) ) ;
    Qi  = blkdiag( [ 0, 1 ; 1, 0 ], -eye( k - 2 ) );
else
    e1i = [] ;
    Ti  = [] ;
    Qi  = [] ;
end

end