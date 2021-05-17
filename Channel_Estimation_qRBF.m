%% START
clc;
close all;
clear all;  

fp = 0.9;
alpha_f=0.5;
alpha_q=0.9;
gamma_q=4;

meu = 1e-2;% Step size
meuq = 1e-2;% Step size
meuf = 1e-2;% Step size
len = 1000; % Length of the signal 
runs = 1000; % Number of times signal passes through ADF for weight adaptation

x=[ones(1,round(len/4)) -ones(1,round(len/4)) ones(1,round(len/4)) -ones(1,round(len/4))];
x=awgn(x,10);
%% Defining Unknown System
h = [2 -0.5 -0.1 -0.7 3];

c = [-5:2:5];
n1=length(c);
q = 1;
q_max = 10;
G=eye(n1);

W = randn(3,n1); % Weights
Wf = W;
Wq = W;

beeta=1;
b=randn(1);
bf = b;
bq = b;
% b=0;

%%
% tic
for k=1:runs
    I(k)=0;
    If(k)=0;
    Iq(k)=0;
    U = zeros(3,1);
    U(2:end)=[-1 -1];
    for i1=1:len
        U(1:end-1)=U(2:end);
        U(end)=x(i1);
        for i2=1:n1
            ED(:,i2)=exp((-(abs(U-c(i2))))/beeta^2);
        end
        
        
        %% RBF
        y(i1)=sum(diag(W*ED'))+b;
        d(i1)= h(1)*U(end) +h(2)*U(end-1)+h(3)*U(end-2)+h(4)*(cos(h(5)*U(end)) +exp(-abs(U(end))))+0.1*randn();
        e=d(i1)-y(i1);
        I(k)=I(k)+e*e'./len;   %%% Objective Function

        W=W+meu*e*ED;
        
        b=b+meu*e;
        
        %% Fractional RBF
        yf(i1)=sum(diag(Wf*ED'))+bf;
        d(i1)= h(1)*U(end) +h(2)*U(end-1)+h(3)*U(end-2)+h(4)*(cos(h(5)*U(end)) +exp(-abs(U(end))))+0.1*randn();
        ef=d(i1)-yf(i1);
        If(k)=If(k)+ef*ef'./len;   %%% Objective Function

        Wf=Wf+meuf*ef*ED.*((1-alpha_f)+alpha_f*abs(Wf.^fp));
        bf=bf+meuf*ef.*((1-alpha_f)+alpha_f*abs(bf.^fp));
        
        
        %% q-RBF
        yq(i1)=sum(diag(Wq*ED'))+bq;
        d(i1)= h(1)*U(end) +h(2)*U(end-1)+h(3)*U(end-2)+h(4)*(cos(h(5)*U(end)) +exp(-abs(U(end))))+0.1*randn();
        eq=d(i1)-yq(i1);
        Iq(k)=Iq(k)+eq*eq'./len;   %%% Objective Function

        Wq=Wq+meuq*eq*ED*q;%*G;
        bq=bq+meuq*eq*q;
        
        q = alpha_q*q + gamma_q*eq^2;
        if (q>q_max)
            q=q_max;
        end
        
        
    end
    q_track(k) = q;
end
% time=toc
save comparison_train.mat

plot(q_track)
figure
semilogy(I,'r')
hold on
semilogy(If,'k')
semilogy(Iq,'b')

figure
len = 200; % Length of the signal 

y=0;
I=0;
d=0;
x=0;
yq=0;
Iq=0;
x=[-1 -1 ones(1,round(len/4)) -ones(1,round(len/4)) ones(1,round(len/4)) -ones(1,round(len/4))];
x=awgn(x,20);

    U(2:end)=[-1 -1];

n1=length(c);
    I=0;
    for i1=1:len-1
        U(1:end-1)=U(2:end);
        U(end)=x(i1);
        for i2=1:n1
            ED(:,i2)=exp((-(abs(U-c(i2))))/beeta^2);
        end
        y(i1)=sum(diag(W*ED'))+b;
        d(i1)= h(1)*U(end) +h(2)*U(end-1)+h(3)*U(end-2)+h(4)*(cos(h(5)*U(end)) +exp(-abs(U(end))));
        e=d(i1)-y(i1);
        SE(i1)=e*e';   %%% Objective Function
        I(i1)=mean(SE);
        
        yf(i1)=sum(diag(Wf*ED'))+bf;
        d(i1)= h(1)*U(end) +h(2)*U(end-1)+h(3)*U(end-2)+h(4)*(cos(h(5)*U(end)) +exp(-abs(U(end))));
        ef=d(i1)-yf(i1);
        SEf(i1)=ef*ef';   %%% Objective Function
        If(i1)=mean(SEf);
        
        yq(i1)=sum(diag(Wq*ED'))+bq;
        d(i1)= h(1)*U(end) +h(2)*U(end-1)+h(3)*U(end-2)+h(4)*(cos(h(5)*U(end)) +exp(-abs(U(end))));
        eq=d(i1)-yq(i1);
        SEq(i1)=eq*eq';   %%% Objective Function
        Iq(i1)=mean(SEq);
    end
    hold on
semilogy(I,'r')
semilogy(If,'k')
semilogy(Iq,'b')
%saveas(gcf,strcat('Iq.png'),'png')
figure

plot(d,'cy')
hold on
plot(y,'r')
plot(yf,'r')
plot(yq,'b')
legend('RBF','FRBF','desired','qRBF')
%saveas(gcf,strcat('RBF.png'),'png')

% save(strcat('MSE','.mat'));
  %save(['q_track','.mat'],'q_track');
%   save(['y','.mat'],'y');   save(['yf','.mat'],'yf'); save(['d','.mat'],'d'); save(['yq','.mat'],'yq');
%  save(['I','.mat'],'I');save(['Iq','.mat'],'Iq');
save comparison_test.mat