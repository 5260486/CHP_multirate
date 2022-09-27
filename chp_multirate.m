%% Parameters
k_c=16.05;      % ��¯ú����λ��ֵ
k_b=3.26;       % ��¯ú����λ��ֵ,MJ/m?

C_b=2600;       %  ��������ϵ���ֹ�¯����ϵ����(MW.S/MPa)

K_h=1;
K1=0.3537;      % ���鷢��Ч�ʣ�30~40%
K_3f=6.328;     % �����¯������ϵ��
K_4f=4.951;

K_f=15;

K21=0.002321;
K22=0.001728;
K3=0.03322;
K4=0.04369;

%% Fuel input 
% ub������ȼ����ֵ
% Qc����¯ú�������¯��������Qb����¯ú�������¯��������Nm?/h����׼������/Сʱ��
u_b1=k_c*Q_c1+k_b*Q_B1; 
u_b2=k_c*Q_c2+k_b*Q_B2;

%% Steam heating balance
all_steam_heat=51.5;
q2=all_steam_heat-q1;
q1=u_b1/u_b2*q2;

%% Overvoltage differential model
% Pd:����ѹ����Pt�����ֻ�ѹ����Qw������ȼ�ϣ�K21��K22������������ϵ��
P_t1=P_d1-K21*(K1*u_b1)^1.3;
P_t2=P_d2-K22*(K1*u_b2)^1.3;

% �������ֻ�����ѹ����
P_t=(K3*P_t1*u_T1+K4*P_t2*u_T2)/(K3*u_T1+K4*u_T2);

%% Boiler heat balance
% ��λʱ���ڹ�¯�����ȱ仯=��λʱ���������¯��������������¯�������� 

function P_d=Pd(K_f,P_t,q_h,u_b)
    P_d=1/C_b*(-K_f*P_t-K_h*q_h+K1*u_b);
end

%% Turbine energy balance
% N_E�����ֻ�������ʣ�MW��Kn�����ֻ�����ϵ����Pt�����ֻ�����ѹ����u_T�������ſ���
function N_E=Ne(N_E,P_t,u_T,Kn)
    N_E=1/K_f*(-N_E+Kn*P_t*u_T);
end
