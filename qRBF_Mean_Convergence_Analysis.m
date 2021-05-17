clc;
clear all;
close all;
% txt='RBFNN_FRBFNN.mat';

epoch=100;%2000;

eta = 1e-3;
eta_f = eta;
alpha = 0.5;
neu = 0.5;
beeta=0.01;

c=[-1:sqrt(beeta):1;-1:sqrt(beeta):1]';

n1 = length(c);
temp=0;
% R=0;


x1=-1:0.1:1;
x2=-1:0.1:1;

% x1=randn(1,1000);
% x2=randn(1,1000);

f=[];
P=[];

for i=1:length(x2)
    f=[f exp(-(x1.^2+x2(i).^2))];
    P=[P; x1' repmat(x2(i),length(x1),1)];
end

[m n] = size(P);

   for i1=1:m
        for i2=1:n1
            phi(i1,i2)=exp((-(norm(P(i1,:)-c(i2,:))^2))/beeta^2);
        end
   end

  
   
% R=phi'*phi;
R= eye(n1);
Rho=f*f';

W_opt = diag(Rho*eye(size(R)))'*inv(R);
W_opt = W_opt./sum(W_opt);
Delta_I = eye(size(R))-eta*(diag(diag(R)));

% R=0;

for run=1:1
    
weight = randn(1,n1);
bias = randn();

weight_f = weight;
bias_f = bias;

    
for k=1:epoch
    I(k)=0;
    I_f(k)=0;
    
%     ind=1:m;%randperm(m);
    ind=randperm(m);
%     k
    for i1=1:m
        for i2=1:n1
            phi(1,i2)=exp((-(norm(P(ind(i1),:)-c(i2,:))^2))/beeta^2);
        end
        
        y(i1)=weight*phi(1,:)';% + bias;
        e=f(ind(i1))-y(i1);
%         I(k)=I(k)+e*e';
               
        weight = weight + eta*e*phi(1,:);
%         bias = bias + eta*e;
        Delta_W =abs(W_opt-weight);
        I(k)=I(k)+mean(Delta_W);

%         temp = 1/(norm(phi(1,:)+1e-10)^2);
%           temp = sum(diag(R));
%         weight = weight + 1*e*phi(i1,:)*temp;
%         bias = bias + 1*e*temp;
        

%         y_f(i1)=weight_f*phi(1,:)' + bias_f;
        
%         e_f=f(ind(i1))-y_f(i1);
%         I_f(k)=I_f(k)+e_f*e_f';
             
%         weight_f = weight_f + eta_f*e_f*phi(i1,:).*(alpha+(1-alpha)*real(weight_f.^(1-neu)));
%         bias_f = bias_f + eta_f*e_f.*(alpha+(1-alpha)*real(bias_f.^(1-neu)));

%         eta_f = 0.9*eta_f + 0.9*e_f^2;
%        if (eta_f>(temp))
%            eta_f = temp;
%        end
       
%         weight_f = weight_f + eta_f*(e_f)*phi(1,:);
%         bias_f = bias_f + eta_f*e_f;
        
%         weight_f = weight_f + eta_f*e_f*phi(1,:)*temp;
%         bias_f = bias_f + eta_f*e_f*temp;

%     R = ((i1-1)*R + phi'*phi)./i1;
    end
    
I_f(k)=mean(abs((Delta_I^(k+1))*(W_opt-weight_f)')) ;      % After approximation 
% I= I./m;
end

 

I100(run,:) = I./m;
I_f100(run,:) = I_f;
end

% R

semilogy(mean(I100,1))
hold on
semilogy(mean(I_f100,1),'r')
% 
% figure
% plot(f,'^-g');
% hold on;
% plot(y(ind),'o-b');
% plot(y_f(ind),'*-r');

