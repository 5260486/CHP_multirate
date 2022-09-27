clc; clear; close all; 

%% Main
format long;
a=0;
b=3600;         %1h
h=1e-3;

%% Parameters
k_c=16.05;      % ��¯ú����λ��ֵ
k_b=3.26;       % ��¯ú����λ��ֵ,MJ/m?

C_b=2600;       %  ��������ϵ���ֹ�¯����ϵ����(MW.S/MPa)

K_h=1;
K1=0.3537;      % ���鷢��Ч�ʣ�30~40%
K_3f=6.328;     % �����¯������ϵ��
K_4f=4.951;

Kf=15;

K21=0.002321;
K22=0.001728;
K3=0.03322;
K4=0.04369;

% ��һ�У�ʱ��Time/h��
% �ڶ��У�Q_B�������У�Q_C��
% �����У�P_b�������У�P_t��
% �����У�u_T�������У�NE��
Boiler1_Data=load('D:\\program\\CHP_experiment\\Boiler1Data.txt');
Boiler2_Data=load('D:\\program\\CHP_experiment\\Boiler2Data.txt');

%% Initial value


%% Fuel input 
% ub������ȼ����ֵ
% Qc����¯ú�������¯��������Qb����¯ú�������¯��������Nm3/h����׼������/Сʱ��
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

    % Pd:����ѹ����Pt�����ֻ�ѹ����Qw������ȼ�ϣ�K21��K22������������ϵ��
    P_t1=P_d1-K21*(K1*u_b1)^1.3;
    P_t2=P_d2-K22*(K1*u_b2)^1.3;

    % �������ֻ�����ѹ����
    P_t=(K3*P_t1*u_T1+K4*P_t2*u_T2)/(K3*u_T1+K4*u_T2);    
    
    % Boiler heat balance
    % ��λʱ���ڹ�¯�����ȱ仯=��λʱ���������¯��������������¯������
    d_P_d1=1/C_b*(-K_3f*P_t-K_h*q1+K1*u_b1);            % d_P_d1=dy1/dt
    d_P_d2=1/C_b*(-K_4f*P_t-K_h*q2+K1*u_b2);
    
end

%% Turbine energy balance
% N_E�����ֻ�������ʣ�MW��Kn�����ֻ�����ϵ����Pt�����ֻ�����ѹ����u_T�������ſ���

function N_nE=Calculate_Ne(N_nE,P_t,u_nT,Kn)        % N_nE=x1,x2
    d_N_nE=1/Kf*(-N_nE+Kn*P_t*u_nT);                      % d_N_nE=dx/dt
end

%% RungeKutta solve

function result=Runge2(a,b,h,P_d_initial)
    k1=Calculate_Pt(K_3f,P_t_initial,q2,u_b1);
    

end

