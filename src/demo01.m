clc
clear
close all

F1 = [1 0; 0 4]; 
P1 = sqrtm(F1);
F2 = [5/2 3/2; 3/2 5/2]; 
P2 = sqrtm(F2);
g1 = [-1; 0]; 
g2 = [-11; -13];
c0 = [0; 0]; 
c1 = P1 \ g1; 
c2 = P2 \ g2;
gamma1 = -3; 
gamma2 = 70; 
d0 = 0;
d1 = sqrt(c1'*c1-gamma1); 
d2 = sqrt(c2'*c2-gamma2);
A0 = [zeros(2,1) eye(2) -eye(2)]';
A1 = [zeros(2,1) P1 zeros(2,2)]';
A2 = [zeros(2,1) zeros(2,2) P2]';
b0 = [1 0 0 0 0]'; 
b1 = zeros(5,1); 
b2 = b1;
At0 = -[b0 A0]; 
At1 = -[b1 A1]; 
At2 = -[b2 A2]
At = [At0 At1 At2];
bt = [-1 0 0 0 0]';
ct0 = [d0; c0]; 
ct1 = [d1; c1]; 
ct2 = [d2; c2]
ct = [ct0; ct1; ct2];
K = []; 
K.q = [size(At0,2) size(At1,2) size(At2,2)];
[xs,ys,info] = sedumi(At,bt,ct,K);
xs
ys
x = ys ;
d = x(1) % Distance: d = norm(u-v,2)
u = x(2:3)
v = x(4:5)

problem.A = At ;
problem.b = bt ;
problem.c = ct ;
problem.K = K ;
save problem problem




