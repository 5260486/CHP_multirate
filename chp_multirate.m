%% Parameters
k_c=16.05;      % 焦炉煤气低位热值
k_b=3.26;       % 高炉煤气低位热值,MJ/m?

C_b=2600;       %  汽包蓄热系数≈锅炉蓄热系数，(MW.S/MPa)

K_h=1;
K1=0.3537;      % 机组发电效率，30~40%
K_3f=6.328;     % 代表锅炉的增益系数
K_4f=4.951;

K_f=15;

K21=0.002321;
K22=0.001728;
K3=0.03322;
K4=0.04369;

%% Fuel input 
% ub：输入燃料热值
% Qc：焦炉煤气输入锅炉的流量；Qb：高炉煤气输入锅炉的流量，Nm?/h（标准立方米/小时）
u_b1=k_c*Q_c1+k_b*Q_B1; 
u_b2=k_c*Q_c2+k_b*Q_B2;

%% Steam heating balance
all_steam_heat=51.5;
q2=all_steam_heat-q1;
q1=u_b1/u_b2*q2;

%% Overvoltage differential model
% Pd:汽包压力；Pt：汽轮机压力；Qw：输入燃料；K21、K22：过热器阻尼系数
P_t1=P_d1-K21*(K1*u_b1)^1.3;
P_t2=P_d2-K22*(K1*u_b2)^1.3;

% 计算汽轮机进气压力：
P_t=(K3*P_t1*u_T1+K4*P_t2*u_T2)/(K3*u_T1+K4*u_T2);

%% Boiler heat balance
% 单位时间内锅炉的蓄热变化=单位时间内流入锅炉的热量与流出锅炉的热量。 

function P_d=Pd(K_f,P_t,q_h,u_b)
    P_d=1/C_b*(-K_f*P_t-K_h*q_h+K1*u_b);
end

%% Turbine energy balance
% N_E：汽轮机输出功率，MW；Kn：汽轮机增益系数；Pt：汽轮机进气压力；u_T进气阀门开度
function N_E=Ne(N_E,P_t,u_T,Kn)
    N_E=1/K_f*(-N_E+Kn*P_t*u_T);
end
