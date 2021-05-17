clc
clear all
close all

load MAE1.mat
I1=mean(I100,1);
If1=mean(I_f100,1);

load MAE2.mat
I2=mean(I100,1);
If2=mean(I_f100,1);

load MAE4.mat
I4=mean(I100,1);
If4=mean(I_f100,1);


figure
plot(I1,':ob','linewidth',2)
hold on
plot(If1,'r','linewidth',2)
plot(I2,':ob','linewidth',2)
plot(If2,'r','linewidth',2)
plot(I4,':ob','linewidth',2)
plot(If4,'r','linewidth',2)
% grid minor
% xlabel('Epoch iterations')
% ylabel('Mean absolute error')
% legend('Simulation','Analytical')
% title('Convergence Performance of q-RBF NN')

leg_handle=legend('Simulation','Analytical');
ax = gca; % current axes
ax.FontSize = 15;
set(leg_handle,'Fontsize',18);
xlabel('Epoch iterations','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Mean Absolute Error','FontSize',16,'FontWeight','bold','Color','k')
grid minor
% ylim([-19 -16])
saveas(gcf,strcat('MAE_Analysis.png'),'png')

corrcoef(If1,I1)
corrcoef(If2,I2)
corrcoef(If4,I4)
