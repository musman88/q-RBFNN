%% Analysis of q-RBF Convergence performance (mean absolute error)
% Author: SHUJAAT KHAN, Shujaat123@gmail.com

clc;
clear all;
close all;

% Mean Convergence Analysis of q-RBF
SNR = 20;
epoch=50;      % Epochs iterations
eta = 0.5e-2;     % Learning rate
beeta=0.1;      % Gaussian kernel spread
neurons = 2;    % Number of neurons
runs = 100;     % Number of independent Monte Carlo simulations
q = 6;
G=q*eye(neurons); % q-Gradient gradient decent-based learning parameter

for run=1:runs

    % Random Input signal
    Input=[beeta*randn(1,100)+0.25 beeta*randn(1,100)+0.75]';
    [m ~] = size(Input);
    
    % K-means clustering to find center of Gaussian kernels
    [~,c,~]=kmeans(Input,neurons);
    
    % Desired output
    output=[zeros(1,100) ones(1,100)]./1000;
    output=awgn(output,SNR);
    
    % Random weights initialization
    weight = randn(1,neurons);
  
    % Simulation (begin)
    for k=1:epoch
    
        I(k)=0; % Initialization of MSE   
        
        % loop to iterate all samples
        for i1=1:m
            
            % Mapping input data to higher dimentions (kernel function)
            for i2=1:neurons
                phi(i1,i2)=exp((-(norm(Input(i1,:)-c(i2,:))^2))/beeta^2);
            end
            
            % Output of the RBF-NN
            y(i1)=weight*phi(i1,:)';
            
            % Instantaneous estimation error
            error=output(i1)-y(i1);
            
            % q-RBF weight update rule
            weight = weight + eta*error*phi(i1,:)*G;
            
            % Accumulating estimation error
            I(k)=I(k)+(error.^2);
    
        end

    end
    
    % Mean absolute error (Simulation) - for each independent run
    I100(run,:) = I./m;
        
end

I = mean(I100,1);
%% Results

% % Simulation Results Log-scale
% figure
% plot(I,'--b','linewidth',2)
% grid minor
% xlabel('Epoch iterations')
% ylabel('Mean absolute error')
% legend('Simulation')
% title('Convergence Performance of q-RBF NN')
% 

save (['I',int2str(q),'.mat'],'I');

Final_Graph_Sensitivity

