clc
clear all
close all

lw = 2;

load('Iq2.mat');
Iq2=Iq;

load('Iq3.mat');
Iq5=Iq;

load('Iq4.mat');
Iq10=Iq;

load('Iq5.mat');
Iq12=Iq;


% load('Iq2.mat');
% Iq2=Iq;
% 
% load('Iq5.mat');
% Iq5=Iq;
% 
% load('Iq10.mat');
% Iq10=Iq;
% 
% load('Iq12.mat');
% Iq12=Iq;

figure
plot(10*log10(Iq2),'mag','linewidth',lw)
hold on
plot(10*log10(Iq5),'b','linewidth',lw)
plot(10*log10(Iq10),'r','linewidth',lw)
plot(10*log10(Iq12),'k','linewidth',lw)
%legend('q-RBF (q=2)','q-RBF (q=5)','q-RBF (q=10)','q-RBF (q=12)');
xlabel('Number of iterations','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Mean Square Error (dB)','FontSize',16,'FontWeight','bold','Color','k')
ax = gca; % current axes
ax.FontSize = 14;
grid minor
ylim([-16 0])
saveas(gcf,strcat('sensitivity.png'),'png')

