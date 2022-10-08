clc; clear; close all; 

a=0;
b=3600;         %1h
h1=50;
h2=5;

n1=(b-a)/h1;
n2=(b-a)/h2;
d=n2/n1;

t0=0;
tspan1=h1;
tspan2=h2;
tol=1e-5;
hmin=0.1;
hmax1=tspan1;
hamx2=tspan2;
%% Parameters
k_c=16.05;      % 焦炉煤气低位热值
k_b=3.26;       % 高炉煤气低位热值,MJ/Nm3

Cal1.C_b=2600;       %  汽包蓄热系数≈锅炉蓄热系数，(MW.S/MPa)

Cal1.K1=0.3537;      % 机组发电效率，30~40%
Cal1.K_3f=6.328;%5.609;    % 代表锅炉的增益系数
Cal1.K_4f=4.951;%4.522;

Cal1.Kf=15;      %单位 s

Cal1.K21=0.002321;%0.001923;%
Cal1.K22=0.001728;%0.001882;%
Cal1.K3=0.03322;%0.02954;%
Cal1.K4=0.04369;%0.04964;%

% 第一列：时刻Time/h，10：00-22：00；
% 第二列：Q_B；第三列：Q_C；
% 第四列：P_d；第五列：P_t；
% 第六列：u_T；第七列：NE；
Boiler1_Data=load('E:\\programs\\matlab\\CHP_multirate-main\\Boiler1Data.txt');
Boiler2_Data=load('E:\\programs\\matlab\\CHP_multirate-main\\Boiler2Data.txt');

%% Initial value
% 19:00, 1 号锅炉和 2 号锅炉高炉煤气与焦炉煤气的流量，以及 1 号汽轮机与 2 号汽轮机的进气阀门开度
time=4;
%time=10;

Q_B1_initial=Boiler1_Data(time,2);            
Q_B2_initial=Boiler2_Data(time,2);
Q_C1_initial=Boiler1_Data(time,3);
Q_C2_initial=Boiler2_Data(time,3);

Q_B.Q_B1=Q_B1_initial;
Q_B.Q_B2=Q_B2_initial;
Q_C.Q_C1=Q_C1_initial;
Q_C.Q_C2=Q_C2_initial;

u_b_initial=Calculat_ub(Q_B,Q_C,k_c,k_b);
q_initial=Calculate_q(u_b_initial);

P_d1_initial=Boiler1_Data(time,4);
P_d2_initial=Boiler2_Data(time,4);

N_E1_initial=Boiler1_Data(time,7);
N_E2_initial=Boiler2_Data(time,7);

% u_T
u_T.u_T1=Boiler1_Data(time,6);
u_T.u_T2=Boiler2_Data(time,6);

%% Variable,20:00 u_b
Q_B.Q_B1=Boiler1_Data(time+1,2);
Q_B.Q_B2=Boiler2_Data(time+1,2);
Q_C.Q_C1=Boiler1_Data(time+1,3);
Q_C.Q_C2=Boiler2_Data(time+1,3);

u_b_variable=Calculat_ub(Q_B,Q_C,k_c,k_b);

%% Main 
Pd1=[];
Pd2=[];
Pt1=[];
Pt2=[];
T=[];
Ne1=[];
Ne2=[];
y1_0=P_d1_initial;
y2_0=P_d2_initial;
x1_0=N_E1_initial;
x2_0=N_E2_initial;
q=q_initial;
u_b=u_b_variable;
tic
for i =1:n1

    result1=RKF45_Pd(tspan1,tol,hmin,hmax1,y1_0,y2_0,u_b,q,u_T,Cal1);           
    y1_0=vpa(result1.y1,5);
    y2_0=vpa(result1.y2,5);
    
    P_t1=vpa(result1.Pt1,5);
    P_t2=vpa(result1.Pt2,5);
    P_t=vpa(result1.Pt,5);
   
    Pd1(end+1)=y1_0;
    Pd2(end+1)=y2_0;
    Pt1(end+1)=P_t1;
    Pt2(end+1)=P_t2;
    T(end+1)=result1.t;
    
    for j=1:d
        result2=Runge2_Ne(h2,x1_0,x2_0,P_t,u_T,Cal1);      %这个步长要小
        x1_0=vpa(result2.Ne1,5);
        x2_0=vpa(result2.Ne2,5);
    end
        Ne1(end+1)=x1_0;
        Ne2(end+1)=x2_0;
end
toc

subplot(3,1,1);
p1=plot(T,Pd1,T,Pd2);grid on;p1(1).LineWidth=2;p1(2).LineWidth=2;
title('13:00--14:00');
% title('19:00--20:00');
ylabel('Pd/Mpa'),legend('Boiler 1','Boiler 2');
subplot(3,1,2);
p2=plot(T,Pt1,T,Pt2);grid on;p2(1).LineWidth=2;p2(2).LineWidth=2;
ylabel('Pt/Mpa'),legend('Turbine 1','Turbine 2');
subplot(3,1,3);
p3=plot(T,Ne1,T,Ne2);grid on;p3(1).LineWidth=2;p3(2).LineWidth=2;
xlabel('Time/h'),ylabel('Ne/MW'),legend('Turbine 1','Turbine 2');

%% Fuel input
% ub：输入燃料能量MJ/h
% Qc：焦炉煤气输入锅炉的流量；Qb：高炉煤气输入锅炉的流量，Nm3/h（标准立方米/小时）

function u_b=Calculat_ub(Q_B,Q_C,k_c,k_b)

    u_b.u_b1=(k_c*Q_C.Q_C1+k_b*Q_B.Q_B1)/3.6/1e3; 
    u_b.u_b2=(k_c*Q_C.Q_C2+k_b*Q_B.Q_B2)/3.6/1e3;

end

%%  Steam heating balance
% q1、q2是两个锅炉产生的热蒸汽量，MW。1MW=3.6*1e3MJ/h
function q=Calculate_q(u_b)
    all_steam_heat=51.5;
    syms q_1 q_2;
    u_b1=u_b.u_b1;
    u_b2=u_b.u_b2;
    
    [q_1,q_2]=solve(q_1+q_2==all_steam_heat,q_1/q_2==u_b1/u_b2);
    q.q1=q_1;
    q.q2=q_2;
end

%% Overvoltage differential model+Boiler heat balance

function [P_t,P_t1,P_t2,d_P_d1,d_P_d2]=Calculate_Pd(P_d1,P_d2,u_b,q,u_T,Cal1)        % P_d1=y1,P_d2=y2

    u_b1=u_b.u_b1;
    u_b2=u_b.u_b2;
    u_T1=u_T.u_T1;
    u_T2=u_T.u_T2;
    
    % Pd:汽包压力,Mpa；Pt：汽轮机压力,MPa；u_b：输入燃料，MJ/Nm3；K21、K22：过热器阻尼系数
    P_t1=P_d1-Cal1.K21*(Cal1.K1*u_b1)^1.3;
    P_t2=P_d2-Cal1.K22*(Cal1.K1*u_b2)^1.3;

    % 计算汽轮机进气压力：
    P_t=(P_t1*(Cal1.K3*u_T1)+P_t2*(Cal1.K4*u_T2))/(Cal1.K3*u_T1+Cal1.K4*u_T2);
    %(Cal1.K4*u_T2*P_t1+Cal1.K3*u_T1*P_t2)/(Cal1.K3*u_T1+Cal1.K4*u_T2);   %(Cal1.K_3f*P_t1+Cal1.K_4f*P_t2)/(Cal1.K_3f+Cal1.K_4f);        

    % Boiler heat balance
    % 单位时间内锅炉的蓄热变化=单位时间内流入锅炉的热量与流出锅炉的热量
    
    d_P_d1=1/Cal1.C_b*(-Cal1.K_3f*P_t1-q.q1+Cal1.K1*u_b1);            % d_P_d1=dy1/dt=f(t,y)     与时间无关
    d_P_d2=1/Cal1.C_b*(-Cal1.K_4f*P_t2-q.q2+Cal1.K1*u_b2);
    
end

%% Turbine energy balance
% N_E：汽轮机输出功率，MW；K3、K4：汽轮机增益系数；Pt：汽轮机进气压力；u_T进气阀门开度,mm

function [d_N_E1,d_N_E2]=Calculate_Ne(N_E1,N_E2,P_t,u_T,Cal1)        % N_nE=x1,x2
    d_N_E1=1/Cal1.Kf*(-N_E1+Cal1.K3*P_t*u_T.u_T1);                      % d_N_nE=dx/dt
    d_N_E2=1/Cal1.Kf*(-N_E2+Cal1.K4*P_t*u_T.u_T2);
end

%% RKF45 solve

function result=RKF45_Pd(tspan,tol,hmin,hmax,y1_0,y2_0,u_b,q,u_T,Cal1)
    flag=true;
    h=hmax;
    t=0;
    
    while flag
        
        [~,~,~,k11,k12]=Calculate_Pd(y1_0,y2_0,u_b,q,u_T,Cal1);
        k11=h*k11;
        k12=h*k12;
        [~,~,~,k21,k22]=Calculate_Pd(y1_0+1/4*k11,y2_0+1/4*k12,u_b,q,u_T,Cal1);
        k21=h*k21;
        k22=h*k22;       
        [~,~,~,k31,k32]=Calculate_Pd(y1_0+3/32*k11+9/32*k21,...
                                        y2_0+3/32*k12+9/32*k22,u_b,q,u_T,Cal1);
        k31=h*k31;
        k32=h*k32;   
        [~,~,~,k41,k42]=Calculate_Pd(y1_0+1932/2197*k11-7200/2197*k21+7296/2197*k31,...
                                        y2_0+1932/2197*k12-7200/2197*k22+7296/2197*k32,u_b,q,u_T,Cal1);
        k41=h*k41;
        k42=h*k42;           
        [~,~,~,k51,k52]=Calculate_Pd(y1_0+439/216*k11-8*k21+3680/513*k31-845*k41/4104,...
                                         y2_0+439/216*k12-8*k22+3680/513*k32-845*k42/4104,...
                                         u_b,q,u_T,Cal1);
        k51=h*k51;
        k52=h*k52;         
        [~,~,~,k61,k62]=Calculate_Pd(y1_0-8/27*k11+2*k21-3544/2565*k31+1859/4104*k41-11/40*k51,...
                                        y2_0-8/27*k12+2*k22-3544/2565*k32+1859/4104*k42-11/40*k52,...
                                        u_b,q,u_T,Cal1);
        k61=h*k61;
        k62=h*k62;  
        % compute |y(i+1)-w(i+1)|/h
        D(1)=k11/360-128/4275*k31-2197/75240*k41+k51/50+2/55*k61;
        D(2)=k12/360-128/4275*k32-2197/75240*k42+k52/50+2/55*k62;
        R=abs(D)/h;
       
        if R<tol            % accurate is acceptable
            t=t+h;
            y1=y1_0+25/216*k11+1408/2565*k31+2197/4104*k41-k51/5+D;
            y2=y2_0+25/216*k12+1408/2565*k32+2197/4104*k42-k52/5+D;
            result.Pd1=y1;
            result.Pd2=y2;
            result.t=t;
            y1_0=y1;
            y2_0=y2;
            [Pt,Pt1,Pt2,~,~]=Calculate_Pd(y1_0,y2_0,u_b,q,u_T,Cal1);
            result.Pt1=Pt1;
            result.Pt2=Pt2;
            result.Pt=Pt;
        end

        % com[ute q
        q=0.84*(tol/R)^(1/4);
        % new step size
        if q<=0.1
            h=0.1*h;
        elseif q>=4
            h=4*h;
        else
            h=h*q;
        end
        
        if h>hmax
            h=hmax;
        end
        
        if t>=tspan
            flag=false;
        elseif t+h>tspan
            h=tspan-t;
            if h<hmin
                flag=false;
            end
        end
        
    end
    
end   

function result=RKF45_Ne(tspan,tol,hmin,hmax,x1_0,x2_0,P_t,u_T,Cal1)


end