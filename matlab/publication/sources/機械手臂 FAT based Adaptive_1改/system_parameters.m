function [m1,m2,l1,l2,lc1,lc2,I1,I2,g,x_wall,k_wall] = system_parameters()
m1 = 0.5; m2 = 0.5; % mass of each linkage
l1 = 0.75; l2 = 0.75; % length of each linkage
lc1 = l1/2; lc2 = l2/2; % center of mass of each linkage
I1 = m1*l1^2/12; I2 = m2*l1^2/12; % moment of inertia of each linkage 
g = 9.81; % gravity
end
