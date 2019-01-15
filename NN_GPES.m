
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% This code is for the Development of Generalized Potenial%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for Carbon Clusters with PBC %%%%%%%%%%%%%%%

% Number of epochs
epochs=100;

% Goal to reach
goal=0;

% Gradient 'cutoff'
min_grad = 1e-10;


mu = 0.001;
mu_dec = 0.1;
mu_inc = 10;
mu_max = 1e10;


data = xlsread('Data_C8_03_16_2010_0K_Training.xlsx');

%no of data pts
Q=length(data);
% Q=1;

%in future clust_size and type will be read from the data file
total_columns=30;%length(data(1,:));

% column 1 in the data file -> cluster size
max_clust_size=8;%max(data(1:Q,1)); %5



% Defining NNs for fc(r), f(r) and f(theta)
S2=25;
net_fc = newff([1.2 4.7],[S2 1],{'logsig' 'purelin'},'trainlm');

% All the 2-body terms
% C-C
load('net_fr_CC_Data_C8_03_16_2010_0K_Training-Trial2.mat');
% net_fr_CC=net2_fr_CC;
% net_fr_CC = newff([1.2 4.7],[S2 1],{'logsig' 'purelin'},'trainlm');
% net_fr_CC.iw{1,1}= [  .400000;  -.400000;   .400000;  -.400135;  -.397665;];
% net_fr_CC.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_fr_CC.b{1,1} = [ -.1880000;  .1821667; -.1763333;  .1704967; .1647197;];
% net_fr_CC.b{2,1}= [-0.1620+0.001]; %#ok<NBRAK>

% net_fr_CC.iw{1,1}= rands(S2,1);%;[  .400000;  -.400000;   .400000;  -.400135;  -.397665;];
% net_fr_CC.lw{2,1}=  rands(1,S2);%[-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_fr_CC.b{1,1}= rands(S2,1);%[ -.1880000;  .1821667; -.1763333;  .1704967; .1647197;];
% net_fr_CC.b{2,1}= rands(1,1);%[-0.1620]; %#ok<NBRAK>


% All the 3-body terms
S3=25;
load('net_ftheta_CCC_Data_C8_03_16_2010_0K_Training-Trial2.mat');
% net_ftheta_CCC=net2_ftheta_CCC;
% net_ftheta_CCC = newff([1.2 4.7; 1.2 4.7 ;1.2 4.7 ],[S3 1],{'logsig' 'purelin'},'trainlm');
% net_ftheta_CCC.iw{1,1}= [.400000  .400000 .400000;  -.400000 -.400000 -.400000;   .400000 .400000 .400000;  -.400135 -.400135 -.400135;  -.397665 -.397665 -.397665;];
% net_ftheta_CCC.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_ftheta_CCC.b{1,1}= [ -.1880000;  .1821667; -.1763333;  .1704967; .1647197;];
% net_ftheta_CCC.b{2,1}= [-0.1620]; %#ok<NBRAK>
% net_ftheta_CCC.iw{1,1}= rands(S3,3);%[.400000  .400000 .400000;  -.400000 -.400000 -.400000;   .400000 .400000 .400000;  -.400135 -.400135 -.400135;  -.397665 -.397665 -.397665;];
% net_ftheta_CCC.lw{2,1}= rands(1,S3);%[-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_ftheta_CCC.b{1,1}= rands(S3,1);%[ -.1880000;  .1821667; -.1763333;  .1704967; .1647197;];
% net_ftheta_CCC.b{2,1}= rands(1,1);%[-0.1620]; %#ok<NBRAK>


S4=5;
net_fphi_CCCC = newff([1.2 4.7; 1.2 4.7 ;1.2 4.7;1.2 4.7 ],[S4 1],{'logsig' 'purelin'},'trainlm');

%%% for scaling purposes only

% minV=min(data(:,total_columns));
% maxV=max(data(:,total_columns));
% delV=maxV-minV;
% 


   

for iQ=1:1:Q
    clust_size(iQ)=data(iQ,1);
    colcount=2;
    for it=1:1:clust_size(iQ)
        for jt=it+1:1:clust_size(iQ)
            
            rs(it,jt,iQ)=data(iQ,colcount);%-1+2*(data(iQ,colcount)-minRFeC)/(maxRFeC-minRFeC); %#ok<AGROW>
            rs(jt,it,iQ)=rs(it,jt,iQ); %#ok<AGROW>
            colcount=colcount+1;
        end
    end
    type(iQ,1:max_clust_size)=[6,6,6,6,6,6,6,6];
    V(iQ)=data(iQ,total_columns);
end


X_fc = getx(net_fc);
X_fr_CC = getx(net_fr_CC);  
X_ftheta_CCC = getx(net_ftheta_CCC); 
X_fphi_CCCC = getx(net_fphi_CCCC);


X_fr=vertcat(X_fr_CC);
X(1:length(X_fr),1)=X_fr;
X_ftheta=vertcat(X_ftheta_CCC);
X(length(X_fr)+1:length(X_fr)+ length(X_ftheta),1)=X_ftheta;
X_fphi=vertcat(X_fphi_CCCC);
X(length(X_fr)+ length(X_ftheta)+1:length(X_fr)+ length(X_ftheta)+ length(X_fphi),1)=X_fphi;





numParameters=length(X_fr)+ length(X_ftheta)+ length(X_fphi);

ii = sparse(1:numParameters,1:numParameters,ones(1,numParameters));%%%

tic
[perf,Ex, Vhat] = calcperf_GPES(net_fc,...
                                net_fr_CC, ...
                                net_ftheta_CCC,...
                                net_fphi_CCCC,...
                                rs, V, Q, clust_size, type);

 toc

tic
% [gXt,jjt,normgX]=calcjejj_GPES(net_fc, ...
%                                 net_fr_CC,...
%                                 net_ftheta_CCC,...
%                                 rs, clust_size,Q,type,Ex);%%%%

toc
time_saved=toc;
for epoch=0:epochs
    tic
    [je,jj,normgX]=calcjejj_GPES(net_fc, ...
                                 net_fr_CC,...
                                 net_ftheta_CCC,...
                                 net_fphi_CCCC,...
                                 rs, clust_size,Q,type,Ex);%%%%
    toc
    normgX;
  
    if(isnan(normgX) == 1)
        normgX;
    end

    % Training Record
    epochPlus1 = epoch+1;
    tr.perf(epochPlus1) = perf;
    tr.mu(epochPlus1) = mu;
    tr.gradient(epochPlus1) = normgX;

    % Stopping Criteria
    %   currentTime = etime(clock,startTime);
    if (perf <= goal)
        stop = 'Performance goal met.';
    elseif (epoch == epochs)
        stop = 'Maximum epoch reached, performance goal was not met.';
        %   elseif (currentTime > time)
        %     stop = 'Maximum time elapsed, performance goal was not met.';
    elseif (normgX < min_grad)
        stop = 'Minimum gradient reached, performance goal was not met.';
    elseif (mu > mu_max)
        stop = 'Maximum MU reached, performance goal was not met.';
        %   elseif (doValidation) & (VV.numFail > max_fail)
        %     stop = 'Validation stop.';
        %   elseif flag_stop
        %     stop = 'User stop.';
    end

    % Progress
    %   if isfinite(show) & (~rem(epoch,show) | length(stop))
    %     fprintf('%s%s%s',this,'-',gradientFcn);
    if isfinite(epochs) fprintf(', Epoch %g/%g',epoch, epochs); end
    %   if isfinite(time) fprintf(', Time %4.1f%%',currentTime/time*100); end
    
    if isfinite(goal) fprintf(', %s %g/%g',upper(net_fc.performFcn),perf,goal); end
    
      
    
    if isfinite(min_grad) fprintf(', Gradient %g/%g',normgX,min_grad); end
    fprintf('\n')
    %   flag_stop=plotperf(tr,goal,this,epoch);
    %     if length(stop) fprintf('%s, %s\n\n',this,stop); end
    %   end

    % Stop when criteria indicate its time
    %   if length(stop)
    %     if (doValidation)
    %     net = VV.net;
    %   end
    %     break
    %   end
    
    % Levenberg Marquardt
    while (mu <= mu_max)
        tic
        % CHECK FOR SINGULAR MATRIX
        [msgstr,msgid] = lastwarn;
        lastwarn('MATLAB:nothing','MATLAB:nothing')
        warnstate = warning('off','all');
        dX = -(jj+ii*mu) \ je;
        [msgstr1,msgid1] = lastwarn;
        flag_inv = isequal(msgid1,'MATLAB:nothing');
        if flag_inv, lastwarn(msgstr,msgid); end;
        warning(warnstate)
        X2 = X + dX;

        net2_fc = net_fc;
        
        X2_fr = X2(1:length(X_fr),1);
        net2_fr_CC = setx(net_fr_CC,X2_fr(1:length(X_fr_CC)));
               
        
        X2_ftheta = X2(length(X_fr)+1:length(X_fr)+ length(X_ftheta),1);
        net2_ftheta_CCC=setx(net_ftheta_CCC,X2_ftheta(1:length(X_ftheta_CCC)));
        
        
        X2_fphi = X2(length(X_fr)+ length(X_ftheta)+1:length(X_fr)+ length(X_ftheta)+length(X_fphi),1);
        net2_fphi_CCCC = setx(net_fphi_CCCC,X2_fphi(1:length(X_fphi_CCCC)));
     

        
        [perf2,Ex, Vhat] = calcperf_GPES(net2_fc, net2_fr_CC, ...
                                                  net2_ftheta_CCC,  ...
                                                  net2_fphi_CCCC,  ...
                                           rs, V, Q, clust_size, type);
%         V-Vhat;
        toc
        tic
        if (perf2 < perf) && flag_inv
%             X = X2; net = net2; %Zb = Zb2; Zi = Zi2; Zl = Zl2;
            
            X = X2; 
            
            net_fc = net2_fc; 
            
            
            net_fr_CC= net2_fr_CC;
                    
            
            net_ftheta_CCC = net2_ftheta_CCC;
            
            
            net_fphiCCCC = net2_fphi_CCCC;
           
        
            %N = N2; Ac = Ac2; El = El2;
            perf = perf2;
            mu = mu * mu_dec;
            if (mu < 1e-20)
                mu = 1e-20;
            end
            break   % Must be after the IF
        end
        mu = mu * mu_inc;
        toc
    end
   
end


