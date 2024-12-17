%  设置第 i 个线性锥的转换矩阵 T 和 Q
% 输入
% k : 锥维数
% 输出
% T : 转换矩阵
% Q : 转换矩阵
function [ e1i, Ti, Qi ] = TransMatLinearCone( k )

if k > 0
    e1i = 1 ;
    Ti  = 1 ;
    Qi  = 1 ;
else
    e1i = [] ;
    Ti  = [] ;
    Qi  = [] ;
end

end

