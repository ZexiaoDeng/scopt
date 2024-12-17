
clc
clear
close all

xi = randi([0, 100], 400, 1) ;

tic
[ Xi, invXi ] = ArrowHeadMat( xi ) ;
toc

tic
pinv( Xi ) ;
toc

norm( ( invXi - pinv( Xi ) ) )



