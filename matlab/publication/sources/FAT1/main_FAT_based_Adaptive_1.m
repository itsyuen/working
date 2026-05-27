% FAT-based force sensorless design without S&L and MRC

clc;clf;clear;

%% initialization
[m1,m2,l1,l2,lc1,lc2,I1,I2,g] = system_parameters();

l = 11; % number of terms of z
T = 17; % period

lambda = diag([10 5]);
Gamma_D = 10*[1*diag([ones(2*l,1)]) diag([zeros(2*l,1)]);1*diag([zeros(2*l,1)]) diag([ones(2*l,1)])];
Gamma_C = 10*[1*diag([ones(2*l,1)]) diag([zeros(2*l,1)]);1*diag([zeros(2*l,1)]) diag([ones(2*l,1)])];
Gamma_g = 5000*diag([ones(2*l,1)]);
Kd = diag([50,200]);

tspan = [0 10]; % simulation time

% initial condition 
X0 = [1 0.35]'; % robot (Cartesian space)
q0 = inverse_kinematics(l1,l2,X0); % transfer into the joint space

W0_reshaped = 1*zeros(1,10*l);
Y0 = [q0(1) q0(2) 0 0 W0_reshaped];


Y=Y0;
%% ode89

[t,Y] = ode89(@(t,Y) FAT_based_Adaptive_1(t,Y,lambda,Gamma_D,Gamma_C,Gamma_g,Kd,l,T),tspan,Y0);
q = [Y(:,1) Y(:,2)]'; % robot trajectory in the joint space
q_dot = [Y(:,3) Y(:,4)]';

%% recall data
n = length(t);
Xd = zeros(2,n); % desired trajectory in the Cartesian space
X = zeros(2,n);  % robot trajectory in the Cartesian space
qd = zeros(2,n); % desired trajectory in the joint space
tau = zeros(2,n); % control input
Dx_hat1 = zeros(1,n);
Dx1 = zeros(1,n);
Dx_hat2 = zeros(1,n);
Dx2 = zeros(1,n);
Dx_hat3 = zeros(1,n);
Dx3 = zeros(1,n);
Dx_hat4 = zeros(1,n);
Dx4 = zeros(1,n);
Cx_hat1 = zeros(1,n);
Cx1 = zeros(1,n);
Cx_hat2 = zeros(1,n);
Cx2 = zeros(1,n);
Cx_hat3 = zeros(1,n);
Cx3 = zeros(1,n);
Cx_hat4 = zeros(1,n);
Cx4 = zeros(1,n);
gx = zeros(2,n);
gx_hat = zeros(2,n);

for i = 1:n
    [Xd(:,i),Xd_dot,Xd_ddot,qd(:,i),qd_dot,qd_ddot] = desired_trajectory_cartesian(t(i),l1,l2);
    [X(:,i),X_dot] = forward_kenimatics(l1,l2,q(:,i),q_dot(:,i));
    [J,J_dot] = Jacobian_matrix(l1,l2,q(:,i),q_dot(:,i));
    [D,C,G,Dx,Cx,Gx] = system_matrix(m1,m2,l1,l2,lc1,lc2,I1,I2,g,q(:,i),q_dot(:,i),J,J_dot);    
    
    e = X(:,i)-Xd(:,i);
    e_dot = X_dot - Xd_dot;
    s = e_dot + lambda*e;
    v = X_dot - lambda*e;
    v_dot = Xd_ddot - lambda*e_dot;

    Z = generate_basis(t(i),l,T);
    WD1_hat = reshape(Y(i,5:5+l-1),[l,1]); % estimation of weighting matrix W
    WD2_hat = reshape(Y(i,5+l:5+2*l-1),[l,1]);
    WD3_hat = reshape(Y(i,5+2*l:5+3*l-1),[l,1]);
    WD4_hat = reshape(Y(i,5+3*l:5+4*l-1),[l,1]);

    WC1_hat = reshape(Y(i,5+4*l:5+5*l-1),[l,1]);
    WC2_hat = reshape(Y(i,5+5*l:5+6*l-1),[l,1]);
    WC3_hat = reshape(Y(i,5+6*l:5+7*l-1),[l,1]);
    WC4_hat = reshape(Y(i,5+7*l:5+8*l-1),[l,1]);

    Wg1_hat = reshape(Y(i,5+8*l:5+9*l-1),[l,1]);
    Wg2_hat = reshape(Y(i,5+9*l:5+10*l-1),[l,1]);
    
    WD_hat = [WD1_hat zeros(l,1); zeros(l,1) WD3_hat ; WD2_hat zeros(l,1) ;zeros(l,1) WD4_hat];
    WC_hat = [WC1_hat zeros(l,1); zeros(l,1) WC3_hat ; WC2_hat zeros(l,1) ;zeros(l,1) WC4_hat];
    Wg_hat = [Wg1_hat zeros(l,1); zeros(l,1) Wg2_hat];

    ZD = [Z zeros(l,1); Z zeros(l,1); zeros(l,1) Z;zeros(l,1) Z];
    ZC = [Z zeros(l,1); Z zeros(l,1); zeros(l,1) Z;zeros(l,1) Z];
    Zg = [Z;Z];
    
    Dx_hat = WD_hat'*ZD;
    Dx1(i) = Dx(1,1);
    Dx_hat1(1,i)=Dx_hat(1,1);
    Dx2(i) = Dx(1,2);
    Dx_hat2(1,i)=Dx_hat(1,2);
    Dx3(i) = Dx(2,1);
    Dx_hat3(1,i)=Dx_hat(2,1);
    Dx4(i) = Dx(2,2);
    Dx_hat4(1,i)=Dx_hat(2,2);
    
    Cx_hat = WC_hat'*ZC;
    Cx1(i) = Cx(1,1);
    Cx_hat1(1,i)=Cx_hat(1,1);
    Cx2(i) = Cx(1,2);
    Cx_hat2(1,i)=Cx_hat(1,2);
    Cx3(i) = Cx(2,1);
    Cx_hat3(1,i)=Cx_hat(2,1);
    Cx4(i) = Cx(2,2);
    Cx_hat4(1,i)=Cx_hat(2,2);

    gx(:,i)=Gx;
    gx_hat(:,i)=Wg_hat'*Zg;

    tau(:,i) = J'*(WD_hat'*ZD*v_dot +WC_hat'*ZC*v +Wg_hat'*Zg -Kd*s);

end


%% plot figures
figure(1);
plot(X(1,1), X(2,1), 'bx','LineWidth',1.5);
hold;
plot(X(1,:), X(2,:),'k','LineWidth',1.5);
plot(Xd(1,:), Xd(2,:),'r--','LineWidth',1.5);
% grid;
axis equal;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Initial Position','Robot','Desired', 'Interpreter','latex','FontSize', 14);
xlabel('x(m)', 'Interpreter','latex','FontSize', 14);
ylabel('z(m)', 'Interpreter','latex','FontSize', 14);
%title('Tracking Result by FAT-based Design', 'Interpreter','latex','FontSize', 24);
ylim([0.34 0.62]);

figure(2);
subplot(2,1,1);
plot(t, q(1,:),'k','LineWidth',1.5);
hold;
plot(t, qd(1,:),'r--','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Actual', 'Desired','Interpreter','latex','FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$q_1$(rad)', 'Interpreter','latex','FontSize', 14);
% title('First joint angle', 'Interpreter','latex','FontSize', 24);
subplot(2,1,2);
plot(t, q(2,:),'k','LineWidth',1.5);
hold;
plot(t, qd(2,:),'r--','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Actual','Desired', 'Interpreter','latex','FontSize', 14,'Location','southeast');
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$q_2$(rad)', 'Interpreter','latex','FontSize', 14);
% title('Second joint angle', 'Interpreter','latex','FontSize', 24);


figure(3);
subplot(2,1,1);
plot(t, tau(1,:),'b','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\tau_1$(N$\cdot$m)', 'Interpreter','latex','FontSize', 14);
%title('Control Input $\tau_1$', 'Interpreter','latex','FontSize', 24);
subplot(2,1,2);
plot(t, tau(2,:),'b','LineWidth',1.5);
% grid;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\tau_2$(N$\cdot$ m)', 'Interpreter','latex','FontSize', 14);
%title('Control Input $\tau_2$', 'Interpreter','latex','FontSize', 24);

%% 
figure(4);
subplot(2,2,1);
plot(t,Dx_hat1,'k','LineWidth',1.5);
hold;
plot(t,Dx1,'r--','LineWidth',1.5);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathbf{D}_{x}(11)$', 'Interpreter','latex','FontSize', 14);
legend('Estimate','Actual', 'Interpreter','latex','FontSize', 20)
subplot(2,2,2);
plot(t,Dx_hat2,'k','LineWidth',1.5);
hold;
plot(t,Dx2,'r--','LineWidth',1.5);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathbf{D}_{x}(12)$', 'Interpreter','latex','FontSize', 14);
legend('Estimate','Actual', 'Interpreter','latex','FontSize', 20)
subplot(2,2,3);
plot(t,Dx_hat3,'k','LineWidth',1.5);
hold;
plot(t,Dx3,'r--','LineWidth',1.5);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathbf{D}_{x}(21)$', 'Interpreter','latex','FontSize', 14);
legend('Estimate','Actual', 'Interpreter','latex','FontSize', 20)
subplot(2,2,4);
plot(t,Dx_hat4,'k','LineWidth',1.5);
hold;
plot(t,Dx4,'r--','LineWidth',1.5);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathbf{D}_{x}(22)$', 'Interpreter','latex','FontSize', 14);
legend('Estimate','Actual', 'Interpreter','latex','FontSize', 20)

%% 

figure(5);
subplot(2,2,1);
plot(t,Cx_hat1,'k','LineWidth',1.5);
hold;
plot(t,Cx1,'r--','LineWidth',1.5);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathbf{C}_{x}(11)$', 'Interpreter','latex','FontSize', 14);
legend('Estimate','Actual', 'Interpreter','latex','FontSize', 20)
subplot(2,2,2);
plot(t,Cx_hat2,'k','LineWidth',1.5);
hold;
plot(t,Cx2,'r--','LineWidth',1.5);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathbf{C}_{x}(12$)', 'Interpreter','latex','FontSize', 14);
legend('Estimate','Actual', 'Interpreter','latex','FontSize', 20)
subplot(2,2,3);
plot(t,Cx_hat3,'k','LineWidth',1.5);
hold;
plot(t,Cx3,'r--','LineWidth',1.5);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathbf{C}_{x}(21)$', 'Interpreter','latex','FontSize', 14);
legend('Estimate','Actual', 'Interpreter','latex','FontSize', 20)
subplot(2,2,4);
plot(t,Cx_hat4,'k','LineWidth',1.5);
hold;
plot(t,Cx4,'r--','LineWidth',1.5);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathbf{C}_{x}(22)$', 'Interpreter','latex','FontSize', 14);
legend('Estimate','Actual', 'Interpreter','latex','FontSize', 20)
%% 

figure(6);
subplot(2,1,1);
plot(t,gx_hat(1,:),'k','LineWidth',1.5);
hold;
plot(t,gx(1,:),'r--','LineWidth',1.5);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathbf{g}_{x}(1)$', 'Interpreter','latex','FontSize', 14);
legend('Estimate','Actual', 'Interpreter','latex','FontSize', 20)
subplot(2,1,2);
plot(t,gx_hat(2,:),'k','LineWidth',1.5);
hold;
plot(t,gx(2,:),'r--','LineWidth',1.5);
xlabel('Time(sec)', 'Interpreter','latex','FontSize', 14);
ylabel('$\mathbf{g}_{x}(2)$', 'Interpreter','latex','FontSize', 14);
legend('Estimate','Actual', 'Interpreter','latex','FontSize', 20)