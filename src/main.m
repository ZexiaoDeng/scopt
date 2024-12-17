clc
clear
close all
format long
load('problem.mat')

A = problem.A ;
b = problem.b ;
c = problem.c ;
K = problem.K ;

pars.eps = 1e-9 ;
[ xs, ys, info ] = sedumi( A, b, c, K, pars ) ;
fprintf( '===================== 我的求解 ===========================\n' ) ;
% [ x, y, s ] = Scopt01( A, b, c, K ) ;
[ x, y, s ] = Scopt02( A, b, c, K ) ;
eps = [ norm( x - xs ), norm( y - ys ) ]
sol = [ y, ys ]