clc; clear; close all; 

%% Main
format long;
a=0;
b=3600;         %1h
h=1e-3;

%% Parameters
k_c=16.05;      % 焦炉煤气低位热值
k_b=3.26;       % 高炉煤气低位热值,MJ/m?

C_b=2600;       %  汽包蓄热系数≈锅炉蓄热系数，(MW.S/MPa)

K_h=1;
K1=0.3537;      % 机组发电效率，30~40%
K_3f=6.328;     % 代表锅炉的增益系数
K_4f=4.951;

Kf=15;

K21=0.002321;
K22=0.001728;
K3=0.03322;
K4=0.04369;

% 第一列：时刻Time/h；
% 第二列：Q_B；第三列：Q_C；
% 第四列：P_b；第五列：P_t；
% 第六列：u_T；第七列：NE；
Boiler1_Data=load('D:\\program\\CHP_experiment\\Boiler1Data.txt');
Boiler2_Data=load('D:\\program\\CHP_experiment\\Boiler2Data.txt');

%% Initial value


%% Fuel input 
% ub：输入燃料热值
% Qc：焦炉煤气输入锅炉的流量；Qb：高炉煤气输入锅炉的流量，Nm3/h（标准立方米/小时）
Q_B1=Boiler1_Data(:,2);
Q_B2=Boiler2_Data(:,2);
Q_C1=Boiler1_Data(:,3);
Q_C2=Boiler2_Data(:,3);

u_b1=k_c*Q_C1+k_b*Q_B1; 
u_b2=k_c*Q_C2+k_b*Q_B2;

%% Steam heating balance
all_steam_heat=51.5;
q1=[];
q2=[];

for i =1:length(Boiler1_Data(:,1))
    syms q_1 q_2;
    [q_1,q_2]=solve(q_1+q_2==all_steam_heat,q_1/q_2==u_b1(i)/u_b2(i));
    q1(end+1)=q_1;
    q2(end+1)=q_2;
end
q1=q1';
q2=q2';

%% Overvoltage differential model+Boiler heat balance

function P_d=Calculate_Pt(P_d1,P_d2)        % P_d1=y1,P_d2=y2

    % Pd:汽包压力；Pt：汽轮机压力；Qw：输入燃料；K21、K22：过热器阻尼系数
    P_t1=P_d1-K21*(K1*u_b1)^1.3;
    P_t2=P_d2-K22*(K1*u_b2)^1.3;

    % 计算汽轮机进气压力：
    P_t=(K3*P_t1*u_T1+K4*P_t2*u_T2)/(K3*u_T1+K4*u_T2);    
    
    % Boiler heat balance
    % 单位时间内锅炉的蓄热变化=单位时间内流入锅炉的热量与流出锅炉的热量
    d_P_d1=1/C_b*(-K_3f*P_t-K_h*q1+K1*u_b1);            % d_P_d1=dy1/dt
    d_P_d2=1/C_b*(-K_4f*P_t-K_h*q2+K1*u_b2);
    
end

%% Turbine energy balance
% N_E：汽轮机输出功率，MW；Kn：汽轮机增益系数；Pt：汽轮机进气压力；u_T进气阀门开度

function N_nE=Calculate_Ne(N_nE,P_t,u_nT,Kn)        % N_nE=x1,x2
    d_N_nE=1/Kf*(-N_nE+Kn*P_t*u_nT);                      % d_N_nE=dx/dt
end

%% RungeKutta solve

function result=Runge2(a,b,h,P_d_initial)
    k1=Calculate_Pt(K_3f,P_t_initial,q2,u_b1);
    

end

