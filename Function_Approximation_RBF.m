clc;
clear all;
close all;
txt='RBFNN_FRBFNN.mat';

x1=-1:.2:1;
x2=-1:.2:1;
f=[];
P=[];

for i=1:length(x2)
    f=[f exp(-(x1.^2+x2(i).^2))];
    P=[P; x1' repmat(x2(i),length(x1),1)];
end

[m n] = size(P);

epoch=100;%2000;

eta = 1e-3;
eta_f = 10*eta;
alpha = 0.5;
neu =0.5;

beeta=1;


c=[0:0.01:1;0:0.01:1]';

n1 = length(c);
temp=0;

for run=1:10
    
weight = randn(1,n1);
bias = 0;%randn();

weight_f = weight;
bias_f = bias;

    
for k=1:epoch
    I(k)=0;
    I_f(k)=0;
    
    ind=randperm(m);
%     k
    for i1=1:m
        for i2=1:n1
            phi(i1,i2)=exp((-(norm(P(ind(i1),:)-c(i2,:))^2))/beeta^2);
        end
        
        y(i1)=weight*phi(i1,:)' + bias;
        e=f(ind(i1))-y(i1);
        I(k)=I(k)+e*e';
        
%         weight = weight + eta*e*phi(i1,:);
%         bias = bias + eta*e;

        temp = 1/(norm(phi(i1,:))^2);
        weight = weight + 1*e*phi(i1,:)*temp;
        bias = bias + 1*e*temp;
        

        y_f(i1)=weight_f*phi(i1,:)' + bias_f;
        
        e_f=f(ind(i1))-y_f(i1);
        I_f(k)=I_f(k)+e_f*e_f';
        
        weight_f = weight_f + eta_f*e_f*phi(i1,:).*(alpha+(1-alpha)*real(weight_f.^(1-neu)));
        bias_f = bias_f + eta_f*e_f.*(alpha+(1-alpha)*real(bias_f.^(1-neu)));

%         weight_f = weight_f + eta_f*e_f*phi(i1,:).*(alpha+(1-alpha)*(abs(weight_f).^(1-neu)));
%         bias_f = bias_f + eta_f*e_f.*(alpha+(1-alpha)*(abs(bias_f).^(1-neu)));
        
%         weight_f = weight_f + eta_f*e_f*phi(i1,:).*(alpha+(1-alpha)*weight_f.*(abs(weight_f).^(-neu)));
%         bias_f = bias_f + eta_f*e_f.*(alpha+(1-alpha)*bias_f.*(abs(bias_f).^(-neu)));


%         eta_f = 0.9*eta_f + 0.9*e_f^2;
%        if (eta_f>(temp))
%            eta_f = temp;
%        end
%        
%         weight_f = weight_f + eta_f*(e_f)*phi(i1,:);
%         bias_f = bias_f + eta_f*e_f;
        
%         weight_f = weight_f + eta_f*e_f*phi(i1,:)*temp;
%         bias_f = bias_f + eta_f*e_f*temp;

    end
    

end

I100(run,:) = I;
I_f100(run,:) = I_f;
end

semilogy(mean(I100,1))
hold on
semilogy(mean(I_f100,1),'r')
% 
% 10*log10([I(end) I_f(end)])
% 
% figure
% plot(f(ind),'^-g');
% hold on;
% plot(y,'o-b');
% plot(y_f,'*-r');




% xt1=-.9:.2:.9;
% xt2=-.9:.2:.9;
% ft=[];
% Pt=[];
% 
% for i=1:length(xt2)
%     ft=[ft exp(-(xt1.^2+xt2(i).^2))];
%     Pt=[Pt; xt1' repmat(xt2(i),1,length(xt1))'];
% end
% 
% [m n] = size(Pt);
% % t=0;
% for i1=1:m
%     for i2=1:n1
% %         t=t+1;
%         phi(i1,i2)=exp((-(norm(Pt(i1,:)-c(i2,:))^2))/beeta^2);
%     end
%     yt(i1)=weight*phi(i1,:)' + bias;
%     yt_f(i1)=weight_f*phi(i1,:)' + bias_f;
% end
% 
% figure
% plot(ft,'^-g');
% hold on;
% plot(yt,'o-b');
% plot(yt_f,'*-r');
% 
% 
% X=[-.9 -.7 -.5 -.3 -.1 .1 .3 .5 .7 .9]';
% Y=[-.9 -.7 -.5 -.3 -.1 .1 .3 .5 .7 .9]';
% 
% for i=1:10
%    Z(i,:)=f(10*i-9:10*i);
%    Zout(i,:)=yt(10*i-9:10*i);
%    Zout_f(i,:)=yt_f(10*i-9:10*i);
% end
% 
% figure
% % subplot(2,1,1)
% % surfc(X,Y,Z)
% % subplot(2,1,2)
% hold on
% % surfc(X,Y,abs(Z-Zout))
% p1=surfc(X,Y,abs(Z-Zout),'FaceColor','b','DisplayName','RBFNN')
% % subplot(2,1,2)
% % hold on
% % surfc(X,Y,abs(Z-Zout_f))
% p2=surfc(X,Y,abs(Z-Zout_f),'FaceColor','r','DisplayName','FRBFNN')
% % colormap(ax1,spring)
% % title('Output for Testing')
% % colorbar
% legend([p1(1) p2(1)],'Absolute error RBF-NN','Absolute error FRBF-NN');
% view([37 12])
% grid minor
% save(txt);
