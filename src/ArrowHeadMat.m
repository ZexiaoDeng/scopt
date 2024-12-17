%  设置第 i 个锥向量的箭头矩阵( Arrow head matrix) mat( x )
% 输入
% xi : 锥向量（LC，QC，RQC）
% 
% 输出
% Xi        : 箭头矩阵
% invXi     : 箭头矩阵的逆阵
% 
function [ Xi, invXi ] = ArrowHeadMat( xi )

n     = length( xi ) ;

if n > 0 & n < 2
    Xi    = 1 ;
    invXi = 1 ;
elseif n >= 2
    
    I     = eye( n - 1 ) ;
    
    SNorm = xi( 2: n )'*xi( 2: n ) ;       % 范数平方

    Xi    = [ xi( 1 )   , xi( 2: n )'         ; ...
              xi( 2: n ), xi(1)*I ; ] ;

    invXi = ( 1/( xi( 1 )^2 - SNorm ) ) ...
          * [  xi( 1 )     , ...
              -xi( 2: n )' ; ...
              -xi( 2: n )  , ...
             ( xi( 1 ) - SNorm/xi( 1 ) )*I + xi( 2: n )*xi( 2: n )'/xi( 1 ) ; ] ;
else
    error( 'Input error for cone vector!' ) ;
end

end




