clear all;
close all;

Ks=80;
Kzeta=diag([1,1]);
Ki=2;
D=0.1*eye(2);
R=0.1*eye(2);
phi=10;

alpha=1.5;
beta=1.5;
p=0.7;
g=10/9;
k=1.3;
v1=1.1;
v2=1/3;
v3=3/4;
v4=1;
rho1=sqrt(2);rho2=sqrt(2);rho3=sqrt(2);rho4=sqrt(2);
d_1=diag([2 2]);

%系统参数
m1=2;
m2=0.85;
L1=0.35;
L2=0.31;
Lc1=0.5;Lc2=0.5;
I1=0.25*m1*L1^2;
I2=0.25*m2*L2^2;
G=9.8;

%状态量
x(1)=pi/2;
x(2)=pi/2;
x(3)=0;
x(4)=0;
dx=[0;0;0;0];

%RBFNN
aij = -5:0.2:5;
c_a =[aij;aij;aij;aij];
Node_a = length(aij);
Sa=zeros(Node_a,1);
Wa=zeros(Node_a,2);
Wa1=5*ones(Node_a,1);
lra=4;
width_a=1;

cij = -5:0.5:5;
Node_c = length(cij);
Sc=zeros(Node_c,1);
Wc=zeros(Node_c,1);
Wc1=5*ones(Node_c,1);
lrc=0.2;
width_c=0.1;

tau1=0;tau2=0;zeta=[0;0];
tau0=[tau1;tau2];
tau=[0;0];
tauL=[-10;-10];
tauH=[10;10];

Time=20;
step_size=0.001;%--------步长
STeps=Time/step_size;%----迭代次数
rec_ratio=1;
rec=2;
outputsize=(STeps-mod(STeps,rec_ratio))/rec_ratio+1;
out(:,1)=[x(1);x(2);x(3);x(4);1;1;0.05;0.1;tau1;tau2;0;0;0;0;0;0;0;0];

for count=2:STeps
    t=(count-1)*step_size;

% 参考轨迹
    x11d=0.1*sin(0.5*t)+cos(0.5*t);
    x12d=(0.1*sin(t)+cos(t));

    x21d=0.05*cos(0.5*t)-0.5*sin(0.5*t);
    x22d=0.1*cos(t)-sin(t);

    ddx11d=-0.025*sin(0.5*t)-0.25*cos(0.5*t);
    ddx12d=-0.1*sin(t)-cos(t);

    x1=[x(1);x(2)];
    x2=[x(3);x(4)];
    xd=[x11d; x12d];
    dxd=[x21d; x22d];
    ddxd=[ddx11d; ddx12d];
    
%     扰动
    d=[0.5*(sin(0.1*t)+1);0.5*(cos(0.1*t)+1)];
    
    e1=x1-xd;    
    e2=x2-dxd;

    %% Critic NN
    for i=1:1:Node_c
        Sc(i,1)=exp(-(e1-cij(:,i))'*(e1-cij(:,i))/width_c^2);
    end
    hatI=Wc'*Sc;
    
    for i=1:1:Node_c
        A(i,1)=-(Sc(i)/phi)+(-2*Sc(i)*(e1-cij(1,i))'/(width_c)^2)*e2;
    end
%     dWc=-lrc*(e1'*D*e1+tau0'*R*tau0+Wc'*A)*A;

    if norm(Wc)<=norm(Wc1)||norm(Wc)==norm(Wc1)&&Wc'*((e1'*D*e1+tau0'*R*tau0+Wc'*A)*A)>0
    dWc=-lrc*(e1'*D*e1+tau0'*R*tau0+Wc'*A)*A;
    elseif norm(Wc)==norm(Wc1)&&Wc'*((e1'*D*e1+tau0'*R*tau0+Wc'*A)*A)<=0
        dWc=-lrc*(e1'*D*e1+tau0'*R*tau0+Wc'*A)*A+lrc*(Wc'*((e1'*D*e1+tau0'*R*tau0+Wc'*A)*A)/(norm(Wc))^2)*Wc;
    end

    
    %% NTSM控制器
    for j=1:2
        ke1a(j)=(alpha*Abs(e1(j),p-inv(k*v1))+beta*Abs(e1(j),g-inv(k*v1)))^(k*v1);
        ke1b(j)=k*v1*(alpha*Abs(e1(j),p-inv(k*v1))+beta*Abs(e1(j),g-inv(k*v1)))^(k*v1-1)*(alpha*(p-inv(k*v1))*Abs(e1(j),p-1/(k*v1))+beta*(g-inv(k*v1))*Abs(e1(j),g-inv(k*v1)));
    end
    Ke1=diag([ke1a(1),ke1a(2)]);
    Ke2=diag([ke1b(1),ke1b(2)]);
    S=Ke1*e1+Sig(e2,v1);
    
    %% Actor NN
    za=[e1;e2];
    for i=1:Node_a
        Sa(i,1)=exp(-(za-c_a(:,i))'*(za-c_a(:,i))/(width_a)^2);
    end
    Fnn=Wa'*Sa;
    prho=Sa*[tanh(Fnn+Ki*hatI)]';
    
    if norm(Wa)<=norm(Wa1)||norm(Wa)==norm(Wa1)&&Wa'*prho>0
        dWa=-lra*prho;
    elseif norm(Wa)==norm(Wa1)&&Wa'*prho<=0
        dWa=-lra*prho+lra*((Wa'*prho)/(norm(Wa))^2)*Wa;
    end
     
    %% System model
    p1=m1*Lc1^2+m2*L1^2+I1;
    p2=m2*Lc2^2+I2;
    p3=m2*L1*Lc2;
    p4=m1*Lc2+m2*L1;
    p5=m2*Lc2;
    
    M1=p1+p2+2*p3*cos(x(2));
    M2=p2+p3*cos(x(2));
    M=[M1 M2;M2 p2];
    M0=[cos(x(2)) cos(x(2));cos(x(2)) 1];

    
    C1=-p3*x(2)*sin(x(2));
    C2=-p3*(x(3)+x(4))*sin(x(2));
    C3=p3*x(3)*sin(x(2));
    C=[C1 C2;C3 0];
    
    G1=p4*G*cos(x(1))+p5*G*cos(x(1)+x(2));
    G2=p5*G*cos(x(1)+x(2));
    G11=[G1;G2];
    
%     L=inv(M)*(-C*x2-G11)+dM*tau;
% 控制输入
    AA1=-inv(v1)*(Ke1+Ke2)*Sig(e2,2-v1);
    BB1=-inv(v1)*diag([Abs(e2(1),1-v1),Abs(e2(2),1-v1)])*(Sig(rho1*Sig(S,v2)+rho2*Sig(S,v3),v4)+zeta+Ks*S);
    tau0=M*(AA1+BB1+ddxd-Fnn-d_1*sign(S));

    %% input saturation
    for i=1:20
    if tau0(1)>=tauH(1)
        tau(1)=tauH(1);
        break;
    elseif tauL(1)<=tau0(1)&&tau0(1)<=tauH(1)
        tau(1)=tau0(1);
        break;
    elseif tau0(1)<=tauL(1)
        tau(1)=tauL(1);
        break;
    end
    end
    
    for i=1:20
   if tau0(2)>=tauH(2)
        tau(2)=tauH(2);
        break;
    elseif tauL(2)<=tau0(2)&&tau0(2)<=tauH(2)
        tau(2)=tau0(2);
        break;
    elseif tau0(2)<=tauL(2)
        tau(2)=tauL(2);
        break;
    end
    end
    dtau=tau-tau0;


    dzeta=-Kzeta*zeta-Sig(rho3*Sig(zeta,v2)+rho4*Sig(zeta,v3),v4)+S+dtau;    %公式（24）
%     dzeta=-Kzeta*zeta-Sig(rho3*Sig(zeta,v2)+rho4*Sig(zeta,v3),v4)+S;

% 状态方程形式
    dx1=x2;
    dx2=inv(M)*(tau0-C*x2-G11+d);
    dx(1)=dx1(1);
    dx(2)=dx1(2);
    dx(3)=dx2(1);
    dx(4)=dx2(2);
    for i=1:4
     x(i)=x(i)+dx(i)*step_size; 
    end

    Wc=Wc+dWc*step_size;
    Wa=Wa+dWa*step_size;
    zeta=zeta+dzeta*step_size;
    
    if S>=0.01
        TIME=t;
    end
    
    if mod(count-1,rec_ratio)==0  
%实际轨迹
        out(1,rec)=x1(1);%x11    关节1位置position  rad
        out(2,rec)=x1(2);%x12     关节2位置position  rad
        out(3,rec)=x2(1);%x21    关节1速度velocity  rad/s
        out(4,rec)=x2(2);%x22     关节2速度velocity  rad/s
% 参考轨迹
        out(5,rec)=x11d;
        out(6,rec)=x12d;
        out(7,rec)=x21d;
        out(8,rec)=x22d;
%不带补偿的输入
        out(9,rec)=tau0(1);
        out(10,rec)=tau0(2);
%输出误差
        out(11,rec)=e1(1);      %位置误差
        out(12,rec)=e1(2);
        out(13,rec)=e2(1);      %速度误差
        out(14,rec)=e2(2);
%输入
        out(15,rec)=tau(1);    %输入力矩1
        out(16,rec)=tau(2);     %输入力矩2
%滑模向量
        out(17,rec)=S(1);
        out(18,rec)=S(2);
 %神经网络辨识
        out(19,rec)=Fnn(1);
        out(20,rec)=Fnn(2);
        rec=rec+1;
    end
end

tout=linspace(0,Time,outputsize-1);
% === 结构赋值并重命名变量 ===
data.tspan    = tout;
data.e_q      = [out(11,:);out(12,:)].';
data.e_q(1,:) = [];
data.e_dq     = [out(13,:);out(14,:)].';
data.tau_mat  = [out(15,:);out(16,:)].';
data.qd_mat   = [out(5,:);out(6,:)].';
data.q_use    = [out(1,:);out(2,:)].';
data.dqd_mat  = [out(7,:);out(8,:)].';
data.dq_use   = [out(3,:);out(4,:)].';
data.rho1   = [out(3,:)];
data.rho2   = [out(3,:)];

% === 保存为结构展开的 .mat 文件（字段名变为变量名）===
save('method2.mat', '-struct', 'data');
% 滑模Sliding vector （m)
figure(1);
plot(tout,out(17,:),'r',tout,out(18,:),'b-.','Linewidth',1.5)

% 实际轨迹与参考轨迹对比（位置）
figure(2);
subplot(2,1,1)
plot(tout,out(1,:),'r',tout,out(5,:),'b-.','Linewidth',1.5);
xlabel('Time[s]','FontSize',14);
ylabel('Angel[rad]','FontSize',14);
LEG=legend('x11','x11_{d}');
set(LEG,'FontName','Times New Roman','FontSize',14);set(gca,'box','off')
subplot(2,1,2)
plot(tout,out(2,:),'r',tout,out(6,:),'b-.','Linewidth',1.5);
xlabel('Time[s]','FontSize',14);
ylabel('Angel[rad]','FontSize',14);
LEG=legend('x12','x12_{d}');
set(LEG,'FontName','Times New Roman','FontSize',14);set(gca,'box','off')

% 实际轨迹与参考轨迹对比（位置）
figure(3);
subplot(2,1,1)
plot(tout,out(3,:),'r',tout,out(7,:),'b-.','Linewidth',1.5);
xlabel('Time[s]','FontSize',14);
ylabel('Angel Speed[rad/s]','FontSize',14);
LEG=legend('x21','x21_{d}');
set(LEG,'FontName','Times New Roman','FontSize',14);set(gca,'box','off')
subplot(2,1,2)
plot(tout,out(4,:),'r',tout,out(8,:),'b-.','Linewidth',1.5);
xlabel('Time[s]','FontSize',14);
ylabel('Angel Speed[rad/s]','FontSize',14);
LEG=legend('x22','x22_{d}');
set(LEG,'FontName','Times New Roman','FontSize',14);set(gca,'box','off')

% 输入力矩
figure(4);
subplot(2,1,1)
plot(tout,out(9,:),'r',tout,out(15,:),'b-.','Linewidth',1.5);
xlabel('Time[s]','FontSize',14);
ylabel('\tau[Voltage]','FontSize',14);
LEG=legend('\tau_0_1','\tau_1');
set(LEG,'FontName','Times New Roman','FontSize',14);set(gca,'box','off')
subplot(2,1,2)
plot(tout,out(10,:),'r',tout,out(16,:),'b-.','Linewidth',1.5);
xlabel('Time[s]','FontSize',14);
ylabel('\tau[Voltage]','FontSize',14);
LEG=legend('\tau_0_2','\tau_2');
set(LEG,'FontName','Times New Roman','FontSize',14);set(gca,'box','off')

% 神经网络辨识结果
figure(5);
plot(tout,out(19,:),'r',tout,out(20,:),'b-.','Linewidth',1.5);
xlabel('Time[s]','FontSize',14);
LEG=legend('Fnn_1','Fnn_2');
set(LEG,'FontName','Times New Roman','FontSize',14);set(gca,'box','off')

% 位置误差
figure(6);
plot(tout,out(11,:),'r',tout,out(12,:),'b','Linewidth',1.5);
xlabel('Time[s]','FontSize',14);
ylabel('rad[Degree]','FontSize',14);
LEG=legend('e_1_1','e_1_2');
set(LEG,'FontName','Times New Roman','FontSize',14);set(gca,'box','off')
axis([0 inf,-1.5 1.5])

% 速度误差
figure(7);
plot(tout,out(13,:),'r',tout,out(14,:),'b-.','Linewidth',1.5);
xlabel('Time[s]','FontSize',14);
ylabel('Position[Degree]','FontSize',14);
LEG=legend('e_2_1','e_2_2');
set(LEG,'FontName','Times New Roman','FontSize',14);set(gca,'box','off')
axis([0 inf,-1.5 1.5])

function y=Sig(x,v)
    y=abs(x).^v.*sign(x);
end

function y1=Abs(x,v)
    y1=abs(x)^v;
end
