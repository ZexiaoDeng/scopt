clc
clear
close all

A = [ 1.5, 3, -0.8, 4 ;
      2, 0, 9, 10 ;
      -7, 4.8, -0.6, 1 ;
      14, 12, 12.3, -4.5 ; ] ;
b = [ 4, 0, 1, -2 ]' ;
[ L, U ] = lu( A ) ;
X = U\( L\ b )

X1 = inv( A )*b

norm( X - X1 )







