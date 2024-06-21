clear;close all;clc;
n=40:40:240;
z=[11e6 5.5e6 2e6 1e6]; 
l_data=8000;
l_head=416;
l_ack=304;

w=32;
m=5; 
L=6;

DIFS=50e-6;
SIFS=10e-6;
sigma=20e-6;
delta=1e-6;
rts=352;
cts=304;

S=[];
Sf=[];
Si=[];
s=[];
s1=[];
s2=[];
Sn=[];
Sn1=[];
St=[];
St1=[];

%saturation throughput for coopmac NOMA based multi-rate WLAN

for i=n
    
    fun=@(x)trans_prob2(x,w,m,L,i);
    x0=[0.2 0.6];
    y=fsolve(fun,x0);
    tau=y(1);

    p_tr=1-(1-tau)^(i/4);
    p_s=((i/4)*tau*(1-tau)^((i/4)-1));
    p_c=p_tr-p_s;
    p_i=1-p_tr 

    t_s=DIFS+((352+304+304)/1e6)+((8184+416)/z(3))+(3*SIFS)+(4*delta);
    t_c=DIFS+(352/1e6)+delta;

    denom=(p_i*sigma)+(p_s*t_s)+(p_c*t_c);
    s2=4*(((p_s*l_data)/denom)/1e6);
    Sn1=[Sn1,s2];

end

plot(n,Sn1,'r-.o');
hold on;
%saturation throughput of NOMA based multi-rate WLAN

for i=n
    
    fun=@(x)trans_prob2(x,w,m,L,i);
    x0=[0.2 0.6];
    y=fsolve(fun,x0);
    tau=y(1);

    p_tr=1-(1-tau)^(i/4);
    p_s=((i/4)*tau*(1-tau)^((i/4)-1));
    p_c=p_tr-p_s;
    p_i=1-p_tr 

    t_s=DIFS+((352+304+304)/1e6)+((8184+416)/z(4))+(3*SIFS)+(4*delta);
    t_c=DIFS+(352/1e6)+delta;

    denom=(p_i*sigma)+(p_s*t_s)+(p_c*t_c);
    s=4*(((p_s*l_data)/denom)/1e6);
    Sn=[Sn,s];

end

plot(n,Sn,'b-.o');
hold on;
%plot(n,Sf,'b--o');
%hold on;

%saturation throughput for coopmac simple multi-rate WLAN

for i=n
    
    fun1=@(x)trans_prob1(x,w,m,L,i);
    x0=[0.1 0.1];
    y=fsolve(fun,x0);
    tau=y(1);
    
    p_tr=1-(1-tau)^i;
    p_s=(i*tau*(1-tau)^(i-1));      
    p_c=p_tr-p_s;
    p_i=1-p_tr;
    t_s1=[];
    t_c1=[];
    
  for j=z
        if j==1e6
            t_s1(4)=(DIFS+((8184+416)/2.75e6)+((304+2*304+352)/1e6)+(5*SIFS)+(4*delta));
            t_c1(4)=(DIFS+((352)/1e6)+delta);
        else
        t_s1=[t_s1,(DIFS+((8184+416)/j)+((304+304+352)/1e6)+(3*SIFS)+(4*delta))];
        t_c1=[t_c1,(DIFS+((352)/1e6)+delta)];
        end
        
    end
    
    idle=p_i*sigma;
    suc=p_s*t_s1(1)+p_s*t_s1(2)+p_s*t_s1(3)+p_s*t_s1(4);
    col=p_c*((1-p_tr)^3)*(t_c1(1)+t_c1(2)+t_c1(3)+t_c1(4));
    col2=(p_tr^2)*((1-p_tr)^2)*(t_c1(2)+t_c1(3)+t_c1(4)+t_c1(3)+t_c1(4)+t_c1(4));
    col3=(p_tr^3)*(1-p_tr)*(t_c1(3)+t_c1(4)+t_c1(4)+t_c1(4)+t_c1(4));
    col4=(p_tr^4)*(max(t_c1));
    
    denom=idle+suc+col+col2+col3+col4;
    S1=4*(((p_s*l_data)/denom)/1e6);
    St1=[St1,S1];
    disp(St1);
    
end
plot(n,St1,'g-.o');
hold on;
%saturation throughput of simple multi-rate WLAN

for i=n
    
    fun1=@(x)trans_prob1(x,w,m,L,i);
    x0=[0.1 0.1];
    y=fsolve(fun,x0);
    tau=y(1);
    
    p_tr=1-(1-tau)^i;
    p_s=(i*tau*(1-tau)^(i-1));      
    p_c=p_tr-p_s;
    p_i=1-p_tr;
    t_s=[];
    t_c=[];
    
    for j=z
        t_s=[t_s,(DIFS+((8184+416)/j)+((304+304+352)/1e6)+(3*SIFS)+(4*delta))];
        t_c=[t_c,(DIFS+((352)/1e6)+delta)];
        
    end
    
    idle=p_i*sigma;
    suc=p_s*t_s(1)+p_s*t_s(2)+p_s*t_s(3)+p_s*t_s(4);
    col=p_c*((1-p_tr)^3)*(t_c(1)+t_c(2)+t_c(3)+t_c(4));
    col2=(p_tr^2)*((1-p_tr)^2)*(t_c(2)+t_c(3)+t_c(4)+t_c(3)+t_c(4)+t_c(4));
    col3=(p_tr^3)*(1-p_tr)*(t_c(3)+t_c(4)+t_c(4)+t_c(4)+t_c(4));
    col4=(p_tr^4)*(max(t_c));
    
    denom=idle+suc+col+col2+col3+col4;
    S=4*(((p_s*l_data)/denom)/1e6);
    St=[St,S];
    disp(St);
    
end

plot(n,St,'k-.o');


legend('Coopmac NOMA WLAN','NOMA based WLAN','Coopmac Conventional WLAN','Conventional WLAN');
ylim([0 10]);
xlim([40 250]);
xlabel('No of users');
ylabel('Aggregate Throughput (Mbps)');
title('Aggregate Throughput vs No of users');