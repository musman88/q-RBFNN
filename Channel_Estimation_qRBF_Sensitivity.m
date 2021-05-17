%% START
clc;
close all;
clear all;  
meu = 1e-2;% Step size
len = 500; % Length of the signal 
epoch = 1000; % Number of times signal passes through ADF for weight adaptation

q=3;

x=[ones(1,round(len/4)) -ones(1,round(len/4)) ones(1,round(len/4)) -ones(1,round(len/4))];
x=awgn(x,10);
%% Defining Unknow System
h = [2 -0.5 -0.1 -0.7 3];

c = [-5:2:5];
n1=length(c);
G=q*eye(n1);

Wq = randn(3,n1); % Weights

beeta=1;
bq=randn(1);

for k=1:epoch
    I(k)=0;
    Iq(k)=0;
    U = zeros(3,1);
    U(2:end)=[-1 -1];
    for i1=1:len
        U(1:end-1)=U(2:end);
        U(end)=x(i1);
        for i2=1:n1
            ED(:,i2)=exp((-(abs(U-c(i2))))/beeta^2);
        end
        
        %% q-RBF
        yq(i1)=sum(diag(Wq*ED'))+bq;
        d(i1)= h(1)*U(end) +h(2)*U(end-1)+h(3)*U(end-2)+h(4)*(cos(h(5)*U(end)) +exp(-abs(U(end))))+0.1*randn();
        eq=d(i1)-yq(i1);
        Iq(k)=Iq(k)+eq*eq'./len;   %%% Objective Function

        Wq=Wq+meu*eq*ED*G;
        bq=bq+meu*eq;
        
        
    end
end

semilogy(Iq)
save(['Iq',int2str(q),'.mat'],'Iq');

Graph_Sensitivity

