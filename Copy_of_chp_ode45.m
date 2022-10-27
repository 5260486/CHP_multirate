clc; clear; close all; 

a=0;
b=3600;         %1h
h1=50;
h2=5;

n1=(b-a)/h1;
n2=(b-a)/h2;
d=n2/n1;

%% Parameters
Cal1.k_c=16.05;      % ???®≤????¶À???
Cal1.k_b=3.26;       % ???®≤????¶À???,MJ/Nm3

Cal1.C_b=2600;       %  ????????????????????????(MW.S/MPa)

Cal1.K1=0.3537;      % ???˙Ä??ßπ???30~40%
Cal1.K_3f=6.328;%5.609;    % ????????????????
Cal1.K_4f=4.951;%4.522;

Cal1.Kf=15;      %??¶À s

Cal1.K21=0.002321;%0.001923;%
Cal1.K22=0.001728;%0.001882;%
Cal1.K3=0.03322;%0.02954;%
Cal1.K4=0.04369;%0.04964;%

% ????ßµ????Time/h??10??00-22??00??
% ????ßµ?Q_B???????ßµ?Q_C??
% ?????ßµ?P_d???????ßµ?P_t??
% ?????ßµ?u_T???????ßµ?NE??
Boiler1_Data=load('D:\\program\\CHP_experiment\\Boiler1Data.txt');
Boiler2_Data=load('D:\\program\\CHP_experiment\\Boiler2Data.txt');

%% Initial value
% 19:00, 1 ?????? 2 ???????®≤?????®≤????????????? 1 ????????? 2 ???????????????????
time=4;
% time=10;

Q_B1_initial=Boiler1_Data(time,2);            
Q_B2_initial=Boiler2_Data(time,2);
Q_C1_initial=Boiler1_Data(time,3);
Q_C2_initial=Boiler2_Data(time,3);

Q_B.Q_B1=Q_B1_initial;
Q_B.Q_B2=Q_B2_initial;
Q_C.Q_C1=Q_C1_initial;
Q_C.Q_C2=Q_C2_initial;

u_b_initial=Calculat_ub(Q_B,Q_C,Cal1);
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

u_b_variable=Calculat_ub(Q_B,Q_C,Cal1);

%% Main 
Pd1=[];
Pd2=[];
Pt1=[];
Pt2=[];
time1=[];
time2=[];
Ne1=[];
Ne2=[];

y1_0=P_d1_initial;
y2_0=P_d2_initial;
x1_0=N_E1_initial;
x2_0=N_E2_initial;
q=q_initial;
u_b=u_b_variable;

m1=0;
m2=0;
p1=h1;
p2=h2;
tspan1=[m1 p1];
tspan2=[m2 p2];

tic
options=odeset('RelTol',1e-5,'AbsTol',1e-7);
for i =1:n1
    [t1,y]=ode15s(@(t,y) Calculate_Pd(t,y,Cal1,u_b,u_T,q),tspan1,[y1_0 y2_0 ]);
    P_d1=y(end,1);
    P_d2=y(end,2);
    P_t1=P_d1-Cal1.K21*(Cal1.K1*u_b.u_b1)^1.3;
    P_t2=P_d2-Cal1.K22*(Cal1.K1*u_b.u_b2)^1.3;
    P_t=(P_t1*(Cal1.K3*u_T.u_T1)+P_t2*(Cal1.K4*u_T.u_T2))/(Cal1.K3*u_T.u_T1+Cal1.K4*u_T.u_T2);
    
    time1(end+1)=p1;
    Pd1(end+1)=P_d1;
    Pd2(end+1)=P_d2;
    Pt1(end+1)=P_t1;
    Pt2(end+1)=P_t2;
    
    m1=m1+h1;
    p1=p1+h1;
    tspan1=[m1 p1];
    y1_0=y(end,1);
    y2_0=y(end,2);
    
    for j=1:d
        [t2,x]=ode15s(@(t,x) Calculate_Ne(t,x,P_t,u_T,Cal1),tspan2,[x1_0 x2_0]);

        time2(end+1)=p2;
        Ne1(end+1)=x(end,1);
        Ne2(end+1)=x(end,2);
        m2=m2+h2;
        p2=p2+h2;
        tspan2=[m2 p2];
        x1_0=x(end,1);
        x2_0=x(end,2);
    end
end
toc
    
subplot(3,1,1);
p1=plot(time1,Pd1,time1,Pd2);grid on;p1(1).LineWidth=2;p1(2).LineWidth=2;
title('13:00--14:00');
% title('19:00--20:00');
ylabel('Pd/Mpa'),legend('Boiler 1','Boiler 2');
subplot(3,1,2);
p2=plot(time1,Pt1,time1,Pt2);grid on;p2(1).LineWidth=2;p2(2).LineWidth=2;
ylabel('Pt/Mpa'),legend('Turbine 1','Turbine 2');
subplot(3,1,3);
p3=plot(time2,Ne1,time2,Ne2);grid on;p3(1).LineWidth=2;p3(2).LineWidth=2;
xlabel('Time/h'),ylabel('Ne/MW'),legend('Turbine 1','Turbine 2');


%% Fuel input
% ub?????????????MJ/h
% Qc?????®≤????????????????Qb?????®≤????????????????Nm3/h???????????/ß≥???

function u_b=Calculat_ub(Q_B,Q_C,Cal1)

    u_b.u_b1=(Cal1.k_c*Q_C.Q_C1+Cal1.k_b*Q_B.Q_B1)/3.6/1e3; 
    u_b.u_b2=(Cal1.k_c*Q_C.Q_C2+Cal1.k_b*Q_B.Q_B2)/3.6/1e3;

end

%%  Steam heating balance
% q1??q2?????????????????????????MW??1MW=3.6*1e3MJ/h
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

function dydt=Calculate_Pd(t,y,Cal1,u_b,u_T,q)

    u_b1=u_b.u_b1;
    u_b2=u_b.u_b2;
    u_T1=u_T.u_T1;
    u_T2=u_T.u_T2;
    
    P_d1=y(1);
    P_d2=y(2);
    
    % Pd:???????,Mpa??Pt??????????,MPa??u_b??????????MJ/Nm3??K21??K22???????????????
    P_t1=P_d1-Cal1.K21*(Cal1.K1*u_b1)^1.3;
    P_t2=P_d2-Cal1.K22*(Cal1.K1*u_b2)^1.3;

    % ??????????????????
    P_t=(P_t1*(Cal1.K3*u_T1)+P_t2*(Cal1.K4*u_T2))/(Cal1.K3*u_T1+Cal1.K4*u_T2);        

    % Boiler heat balance
    % ??¶À?????????????Å£=??¶À????????????????????????????????
    
    d_P_d1=1/Cal1.C_b*(-Cal1.K_3f*P_t1-q.q1+Cal1.K1*u_b1);            % d_P_d1=dy1/dt=f(t,y)     ????????
    d_P_d2=1/Cal1.C_b*(-Cal1.K_4f*P_t2-q.q2+Cal1.K1*u_b2);  

    dydt=double([d_P_d1;d_P_d2]);
end

%% Turbine energy balance
% N_E???????????????MW??K3??K4????????????????Pt????????????????u_T???????????,mm

function dxdt=Calculate_Ne(t,x,P_t,u_T,Cal1)        % N_nE=x1,x2
    N_E1=x(1);
    N_E2=x(2);   
    d_N_E1=1/Cal1.Kf*(-N_E1+Cal1.K3*P_t*u_T.u_T1);                      % d_N_nE=dx/dt
    d_N_E2=1/Cal1.Kf*(-N_E2+Cal1.K4*P_t*u_T.u_T2);
    dxdt=double([d_N_E1;d_N_E2]);
end

%% RKF45 solve

% function result=RKF45_Pd(tspan,tol,hmin,hmax,...
%                                             u_b,q,u_T,Cal1,Ay,...
%                                             y1_0,y2_0,x1_0,x2_0)
%     t=0;
%     h=hmax;
%     flag=true;
%     
%     Pd1=[];
%     Pd2=[];
%     Pt1=[];
%     Pt2=[];
%     Ne1=[];
%     Ne2=[];
%     T=[];
%     
%     while flag
%         
%         result_pd=Calculate_Pd(y1_0,y2_0,x1_0,x2_0,u_b,q,u_T,Cal1);
%         k11=h*result_pd.dpd1;
%         k12=h*result_pd.dpd2;
%         k13=h*result_pd.dNe1;
%         k14=h*result_pd.dNe2;
%         result_pd=Calculate_Pd(y1_0+Ay.Ay1(1)*k11,y2_0+Ay.Ay1(1)*k12,...
%                                        x1_0+Ay.Ay1(1)*k13,x2_0+Ay.Ay1(1)*k14,...
%                                        u_b,q,u_T,Cal1);
%         k21=h*result_pd.dpd1;
%         k22=h*result_pd.dpd2;
%         k23=h*result_pd.dNe1;
%         k24=h*result_pd.dNe2;      
%         result_pd=Calculate_Pd(y1_0+Ay.Ay2(1)*k11+Ay.Ay2(2)*k21,y2_0+Ay.Ay2(1)*k12+Ay.Ay2(2)*k22,...
%                                         x1_0+Ay.Ay2(1)*k13+Ay.Ay2(2)*k23,x2_0+Ay.Ay2(1)*k14+Ay.Ay2(2)*k24,...
%                                         u_b,q,u_T,Cal1);
%         k31=h*result_pd.dpd1;
%         k32=h*result_pd.dpd2;
%         k33=h*result_pd.dNe1;
%         k34=h*result_pd.dNe2;   
%         result_pd=Calculate_Pd(y1_0+Ay.Ay3(1)*k11+Ay.Ay3(2)*k21+Ay.Ay3(3)*k31,...
%                                         y2_0+Ay.Ay3(1)*k12+Ay.Ay3(2)*k22+Ay.Ay3(3)*k32,...
%                                         x1_0+Ay.Ay3(1)*k13+Ay.Ay3(2)*k23+Ay.Ay3(3)*k33,...
%                                         x2_0+Ay.Ay3(1)*k14+Ay.Ay3(2)*k24+Ay.Ay3(3)*k34,...
%                                         u_b,q,u_T,Cal1);
%         k41=h*result_pd.dpd1;
%         k42=h*result_pd.dpd2;
%         k43=h*result_pd.dNe1;
%         k44=h*result_pd.dNe2;            
%         result_pd=Calculate_Pd(y1_0+Ay.Ay4(1)*k11+Ay.Ay4(2)*k21+Ay.Ay4(3)*k31+Ay.Ay4(4)*k41,...
%                                         y2_0+Ay.Ay4(1)*k12+Ay.Ay4(2)*k22+Ay.Ay4(3)*k32+Ay.Ay4(4)*k42,...
%                                         x1_0+Ay.Ay4(1)*k13+Ay.Ay4(2)*k23+Ay.Ay4(3)*k33+Ay.Ay4(4)*k43,...
%                                         x2_0+Ay.Ay4(1)*k14+Ay.Ay4(2)*k24+Ay.Ay4(3)*k34+Ay.Ay4(4)*k44,...
%                                         u_b,q,u_T,Cal1);
%         k51=h*result_pd.dpd1;
%         k52=h*result_pd.dpd2;
%         k53=h*result_pd.dNe1;
%         k54=h*result_pd.dNe2;        
%         result_pd=Calculate_Pd(y1_0+Ay.Ay5(1)*k11+Ay.Ay5(2)*k21+Ay.Ay5(3)*k31+Ay.Ay5(4)*k41+Ay.Ay5(5)*k51,...
%                                         y2_0+Ay.Ay5(1)*k12+Ay.Ay5(2)*k22+Ay.Ay5(3)*k32+Ay.Ay5(4)*k42+Ay.Ay5(5)*k52,...
%                                         x1_0+Ay.Ay5(1)*k13+Ay.Ay5(2)*k23+Ay.Ay5(3)*k33+Ay.Ay5(4)*k43+Ay.Ay5(5)*k53,...
%                                         x2_0+Ay.Ay5(1)*k14+Ay.Ay5(2)*k24+Ay.Ay5(3)*k34+Ay.Ay5(4)*k44+Ay.Ay5(5)*k54,...
%                                         u_b,q,u_T,Cal1);
%         k61=h*result_pd.dpd1;
%         k62=h*result_pd.dpd2;
%         k63=h*result_pd.dNe1;
%         k64=h*result_pd.dNe2;
%         
%         % compute |y(i+1)-w(i+1)|/h
%         D(1)=k11/360-128/4275*k31-2197/75240*k41+k51/50+2/55*k61;
%         D(2)=k12/360-128/4275*k32-2197/75240*k42+k52/50+2/55*k62;
%         D(3)=k13/360-128/4275*k33-2197/75240*k43+k53/50+2/55*k63;
%         D(4)=k14/360-128/4275*k34-2197/75240*k44+k54/50+2/55*k64;
%         R=vpa(norm(D,2)/h,5);
%        
%         if R<tol            % accurate is acceptable
%             t=t+h;
%             y1=vpa(y1_0+25/216*k11+1408/2565*k31+2197/4104*k41-k51/5+D(1),5);
%             y2=vpa(y2_0+25/216*k12+1408/2565*k32+2197/4104*k42-k52/5+D(2),5);
%             x1=vpa(x1_0+25/216*k13+1408/2565*k33+2197/4104*k43-k53/5+D(1),5);
%             x2=vpa(x2_0+25/216*k14+1408/2565*k34+2197/4104*k44-k54/5+D(2),5);    
%             
%             y1_0=y1;
%             y2_0=y2;
%             result_pd=Calculate_Pd(y1_0,y2_0,x1_0,x2_0,u_b,q,u_T,Cal1);
%             Pt_1=result_pd.Pt1;
%             Pt_2=result_pd.Pt2;
%             
% %             result.Pd1=y1;
% %             result.Pd2=y2;
% %             result.Ne1=x1;
% %             result.Ne2=x2;
% %             result.Pt1=vpa(Pt1,5);
% %             result.Pt2=vpa(Pt2,5);
% %             result.Pt=vpa(Pt,5);
% %             result.t=t;
%             Pt1(end+1)=Pt_1;
%             Pt2(end+1)=Pt_2;
%             Pd1(end+1)=y1;
%             Pd2(end+1)=y2;
%             Ne1(end+1)=x1;
%             Ne2(end+1)=x2;
%             T(end+1)=t;
%                        
%         end
% 
%         % compute e
%         e=0.84*(tol/R)^(1/4);
%         % new step size
%         if e<=0.1
%             h=0.1*h;
%         elseif e>=4
%             h=4*h;
%         else
%             h=h*e;
%         end
%         
%         if h>hmax
%             h=hmax;
%         end
%         
%         if t>=tspan
%             flag=false;
%         elseif t+h>tspan
%             h=tspan-t;
%             if h<hmin
%                 flag=false;
%             end
%         end
%     
%         if flag==false
%             subplot(3,1,1);
%             p1=plot(T,Pd1,T,Pd2);grid on;p1(1).LineWidth=2;p1(2).LineWidth=2;
%             title('13:00--14:00');
%             % title('19:00--20:00');
%             ylabel('Pd/Mpa'),legend('Boiler 1','Boiler 2');
%             subplot(3,1,2);
%             p2=plot(T,Pt1,T,Pt2);grid on;p2(1).LineWidth=2;p2(2).LineWidth=2;
%             ylabel('Pt/Mpa'),legend('Turbine 1','Turbine 2');
%             subplot(3,1,3);
%             p3=plot(T,Ne1,T,Ne2);grid on;p3(1).LineWidth=2;p3(2).LineWidth=2;
%             xlabel('Time/h'),ylabel('Ne/MW'),legend('Turbine 1','Turbine 2');
%         end
%             
%     end
% end   

