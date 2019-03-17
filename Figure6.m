% Olsen, Roepstorff & Bang: "Knowing whom to learn from: individual
% differences in metacognition and weighting of social information"
%
% Dan Bang danbang.db@gmail.com 2019

%% -----------------------------------------------------------------------
%% PREPARATION
%% -----------------------------------------------------------------------

% fresh memory
clc;clear;close;

% add paths
addpath('helpers');

% subjects
n_subjects= 39;

% load data
load('data.mat');

%% -----------------------------------------------------------------------
%% ANALYSIS
%% -----------------------------------------------------------------------

% initial response: 1
% revised response: 2
% bad advisor: B (advisor 1)
% good advisor: G (advisor 2)

% loop through subjects
for i_sbj= 1:n_subjects
    
    % preprocess
    Btmp= [];
    fnames= fieldnames(data);
    for i_f = 1:length(fnames); 
        eval(['Btmp.',fnames{i_f},'=[data.',fnames{i_f},'(data.sbjID==i_sbj & data.advisor==1)];']); 
    end
    Gtmp= [];
    fnames= fieldnames(data);
    for i_f = 1:length(fnames); 
        eval(['Gtmp.',fnames{i_f},'=[data.',fnames{i_f},'(data.sbjID==i_sbj & data.advisor==2)];']); 
    end
    
    % HISTORY EFFECTS ON INITIAL RESPONSE
    % specify indices for (c) current and (p) previous trial
    c_idx= 2:120;
    p_idx= 1:119;
    % specify variables of interest
    % here we do analysis in UNSIGNED space
    % bad partner data
    Bc_con1=zscore(abs(Btmp.sbjCon1(c_idx))); % subject's initial confidence on c trial
    Bc_con2=zscore(abs(Btmp.sbjCon2(c_idx))); % subject's revised confidence on c trial
    Bp_con2=zscore(abs(Btmp.sbjCon2(p_idx))); % subject's revised confidence on p trial
    Bp_pcon=zscore(abs(Btmp.advCon(p_idx)));  % partner's confidence on p trial
    Bc_sacc=Btmp.sbjAcc2(c_idx)-.5; % subject's accuracy (now -.5 and +.5) on c trial 
    Bc_pacc=Btmp.advAcc(c_idx)-.5;  % partner's accuracy (now -.5 and +.5) on c trial 
    Bp_sacc=Btmp.sbjAcc2(p_idx)-.5; % subject's accuracy (now -.5 and +.5) on p trial 
    Bp_pacc=Btmp.advAcc(p_idx)-.5;  % partner's accuracy (now -.5 and +.5) on p trial 
    % good partner data
    Gc_con1=zscore(abs(Gtmp.sbjCon1(c_idx))); % subject's initial confidence on c trial
    Gc_con2=zscore(abs(Gtmp.sbjCon2(c_idx))); % subject's revised confidence on c trial
    Gp_con2=zscore(abs(Gtmp.sbjCon2(p_idx))); % subject's revised confidence on p trial
    Gp_pcon=zscore(abs(Gtmp.advCon(p_idx)));  % partner's confidence on p trial
    Gc_sacc=Gtmp.sbjAcc2(c_idx)-.5; % subject's accuracy (now -.5 and +.5) on c trial 
    Gc_pacc=Gtmp.advAcc(c_idx)-.5;  % partner's accuracy (now -.5 and +.5) on c trial 
    Gp_sacc=Gtmp.sbjAcc2(p_idx)-.5; % subject's accuracy (now -.5 and +.5) on p trial 
    Gp_pacc=Gtmp.advAcc(p_idx)-.5;  % partner's accuracy (now -.5 and +.5) on p trial 
    % specify (Y) variable which we would like to predict
    Y= [abs(Bc_con1) abs(Gc_con1)]; % my initial confidence
    % specify (X) variables which we would like to use as predictors
    X= [Bp_con2 Gp_con2; ... % my previous confidence
        Bp_sacc Gp_sacc; ... % my previous accuracy
        Bp_sacc.*Bp_con2 Gp_sacc.*Gp_con2; ... % my previous confidence x accuracy
        Bp_pcon Gp_pcon; ... % your previous confidence
        Bp_pacc Gp_pacc; ... % your previous accuracy
        Bp_pacc.*Bp_pcon Gp_pacc.*Gp_pcon; ... % your previous confidence x accuracy
        ];    
    % run multiple-linear regression (here ' transposes variables)
    [B,DEV,STATS]= glmfit(X',Y','normal');
    % save beta-weights where k is subject index
    sbeta1(i_sbj,:) = B(2:end); % we are for now not interested in intercept
    
    % HISTORY EFFECTS ON REVISED RESPONSE
    % specify indices for (c) current and (p) previous trial
    c_idx= 2:120;
    p_idx= 1:119;
    % specify variables of interest
    % here we do analysis in UNSIGNED space
    % bad partner data
    Bc_con1=zscore(abs(Btmp.sbjCon1(c_idx))); % subject's initial confidence on c trial
    Bc_con2=zscore(abs(Btmp.sbjCon2(c_idx))); % subject's revised confidence on c trial
    Bp_con2=zscore(abs(Btmp.sbjCon2(p_idx))); % subject's revised confidence on p trial
    Bc_pcon=zscore(abs(Btmp.advCon(c_idx))); % partner's confidence on c trial
    Bp_pcon=zscore(abs(Btmp.advCon(p_idx))); % partner's confidence on p trial
    Bc_sacc=Btmp.sbjAcc1(c_idx)-.5; % partner's accuracy (now -.5 and +.5) on c trial 
    Bc_pacc=Btmp.advAcc(c_idx)-.5; % partner's accuracy (now -.5 and +.5) on c trial 
    Bp_sacc=Btmp.sbjAcc2(p_idx)-.5; % partner's accuracy (now -.5 and +.5) on p trial 
    Bp_pacc=Btmp.advAcc(p_idx)-.5; % partner's accuracy (now -.5 and +.5) on p trial 
    % good partner data
    Gc_con1=zscore(abs(Gtmp.sbjCon1(c_idx))); % subject's initial confidence on c trial
    Gc_con2=zscore(abs(Gtmp.sbjCon2(c_idx))); % subject's revised confidence on c trial
    Gp_con2=zscore(abs(Gtmp.sbjCon2(p_idx))); % subject's revised confidence on p trial
    Gc_pcon=zscore(abs(Gtmp.advCon(c_idx))); % partner's confidence on c trial
    Gp_pcon=zscore(abs(Gtmp.advCon(p_idx))); % partner's confidence on p trial
    Gc_sacc=Gtmp.sbjAcc1(c_idx)-.5; % partner's accuracy (now -.5 and +.5) on c trial 
    Gc_pacc=Gtmp.advAcc(c_idx)-.5; % partner's accuracy (now -.5 and +.5) on c trial 
    Gp_sacc=Gtmp.sbjAcc2(p_idx)-.5; % partner's accuracy (now -.5 and +.5) on p trial 
    Gp_pacc=Gtmp.advAcc(p_idx)-.5; % partner's accuracy (now -.5 and +.5) on p trial 
    % specify (Y) variable which we would like to predict
    Y= [Bc_con2 Gc_con2]; % my revised confidence
    % specify (X) variables which we would like to use as predictors
    X= [Bc_con1 Gc_con1; ... % my initial confidence
        Bp_sacc Gp_sacc; ... % my previous accuracy
        Bp_con2 Bp_con2; ... % my previous confidence
        Bp_con2.*Bp_sacc Gp_con2.*Gp_sacc; % my previous confidenc modulated by my previous accuracy
        Bc_con1.*Bp_con2 Gc_con1.*Gp_con2; ... % my initial confidence modulated by my previous confidence
        Bc_con1.*Bp_sacc Gc_con1.*Gp_sacc; ... % my initial confidence modulated by my previous accuracy
        Bc_con1.*Bp_con2.*Bp_sacc Gc_con1.*Gp_con2.*Gp_sacc; ... % my initial confidence modulated by my previous confidence x accuracy
        Bc_pcon Gc_pcon; ... % your advice
        Bp_pacc Gp_pacc; ... % your previous accuracy
        Bp_pcon Gp_pcon; ... % your previous accuracy
        Bp_pcon.*Bp_pacc Gp_pcon.*Gp_pacc; % your previous confidenc modulated by your previous accuracy
        Bc_pcon.*Bp_pcon Gc_pcon.*Gp_pcon; ... % your advice modulated by your previous confidence
        Bc_pcon.*Bp_pacc Gc_pcon.*Gp_pacc; ... % your advice modulated by your previous accuracy
        Bc_pcon.*Bp_pcon.*Bp_pacc Gc_pcon.*Gp_pcon.*Gp_pacc; ... % your advice modulated by your previous confidence x accuracy
        ];    
    % run multiple-linear regression (here ' transposes variables)
    % save beta-weights where k is subject index
    [B,DEV,STATS]= glmfit(X',Y','normal');
    sbeta2(i_sbj,:) = B(2:end); % we are for now not interested in intercept and my initial confidence
    
end

%% -----------------------------------------------------------------------
%% FIGURE
%% -----------------------------------------------------------------------

% HISTORY EFFECTS ON INITIAL CONFIDENCE
% label the predictors that we used above
% here C: confidence; A: accuracy; t: current trial; s: subject; p: partner
% the weights tell us the influence of the predictors on my revised confidence
predictor_names= {'C_s_,_t_-_1','A_s_,_t_-_1','C_s_,_t_-_1xA_s_,_t_-_1', ...
                  'C_p_,_t_-_1','A_p_,_t_-_1','C_p_,_t_-_1xA_p_,_t_-_1'};
n_x= length(predictor_names); % number of predictors
% figure
figz=figure('color',[1 1 1]);
% colours
me_col= [255,105,180]./255;
you_col= [0,255,127]./255;
% simple statistics for plotting
n_s= size(sbeta1,1); % number of subjects
mean1= mean(sbeta1); % mean beta-weights
sem1= std(sbeta1)./(sqrt(n_s)); % SEM beta-weights
% plot reference line going through zero on y-axis
plot([0 0],[0 n_x+1],'k-','LineWidth',2); hold on;
% plot error bars
o= 0; % offset
i_v=1:n_x; for i=1:length(i_v); plot([-1 +1],[i_v(i)-o i_v(i)-o],'k--','LineWidth',.1); hold on; end;
% me
i_v=1:3; for i=1:length(i_v); plot([mean1(i_v(i))-sem1(i_v(i)) mean1(i_v(i))+sem1(i_v(i))],[i_v(i)-o i_v(i)-o],'k-','LineWidth',2); hold on; end;
i_v=1:3; for i=1:length(i_v); plot(mean1(i_v(i)),i_v(i)-o,'ko','LineWidth',1,'MarkerFaceColor',me_col,'MarkerSize',12,'LineWidth',2); hold on; end;
% you
i_v=4:6; for i=1:length(i_v); plot([mean1(i_v(i))-sem1(i_v(i)) mean1(i_v(i))+sem1(i_v(i))],[i_v(i)-o i_v(i)-o],'k-','LineWidth',2); hold on; end;
i_v=4:6; for i=1:length(i_v); plot(mean1(i_v(i)),i_v(i)-o,'ko','LineWidth',1,'MarkerFaceColor',you_col,'MarkerSize',12,'LineWidth',2); hold on; end;
% general settings
set(gca,'YAxisLocation','left','FontSize',28,'LineWidth',2); % set font size for figure
box off; % only lines on x- and y-axis
% specify x-axis
set(gca,'XTick',[-.05:.05:.15]); % specify ticks (no need for labels as numeric)
xlabel('beta-weight','FontSize',36); % specify x-axis name
xlim([-.05 .15]); % set x-axis limits
% specify y-axis
set(gca,'YTick',1:n_x,'YTickLabel',predictor_names); % specify ticks and tick labels
ylabel('predictor','FontSize',36); % specify x-axis name
ylim([0 n_x+1]); % set y-axis limits
print('-djpeg','-r300',['Plots',filesep,'Figure6A']);
% statistics
[a1,b1,c1,d1]= ttest(sbeta1,0);
% print statistics
fprintf(['=================================================\n']);
fprintf(['HISTORY EFFECTS ON INITIAL CONFIDENE\n']);
fprintf(['=================================================\n']);
fprintf(['p-vals: ' num2str(b1),'\n']);
fprintf(['t-vals: ' num2str(d1.tstat),'\n']);
fprintf(['=================================================\n']);

% VISUALISE HISTORY EFECTS ON INITIAL CONFIDNECE
% calculate expected impact
muBeta= mean(sbeta1);
myConPre= muBeta(1);
myAccPre= muBeta(2);
myConXAccPre= muBeta(3);
c_v= -3:1:3;
a_v= [-.5 .5];
for i_c= 1:length(c_v)
    for i_a= 1:2;
        impact(i_a,i_c)= myConPre*c_v(i_c) + myAccPre*a_v(i_a) + myConXAccPre*c_v(i_c)*a_v(i_a);
    end
end
% figure
figure('color',[1 1 1]);
hold on;
plot(c_v,impact(1,:),'r-','LineWidth',2);
plot(c_v,impact(2,:),'g-','LineWidth',2);
% set axes
set(gca,'FontSize',28,'LineWidth',2);
set(gca,'XTick',c_v);
set(gca,'YTick',-.2:.1:.2);
xlim([-3 3]);
ylim([-.2 .2]);
xlabel('C_s_,_t_-_1 [SD]','FontSize',36);
ylabel('impact','FontSize',36);
print('-djpeg','-r300',['Plots',filesep,'Figure6B']);

% HISTORY EFFECTS ON FINAL CONFIDENCE
% label the predictors that we used above
% here C: confidence; A: accuracy; t: current trial; s: subject; p: partner
% the weights tell us the influence of the predictors on my revised confidence
predictor_names= {'C_s_,_t','A_s_,_t_-_1','C_s_,_t_-_1','A_s_,_t_-_1xC_s_,_t_-_1','C_s_,_txC_s_,_t_-_1','C_s_,_txA_s_,_t_-_1','C_s_,_txC_s_,_t_-_1xA_s_,_t_-_1', ...
                  'C_p_,_t','A_p_,_t_-_1','C_s_,_t_-_1','A_p_,_t_-_1xC_p_,_t_-_1','C_p_,_txC_p_,_t_-_1','C_p_,_txA_p_,_t_-_1','C_p_,_txC_p_,_t_-_1xA_p_,_t_-_1'};
n_x= length(predictor_names); % number of predictors
% figure
figz=figure('color',[1 1 1]);
% colours
me_col= [255,105,180]./255;
you_col= [0,255,127]./255;
% simple statistics for plotting
n_s= size(sbeta2,1); % number of subjects
mean1= mean(sbeta2); % mean beta-weights
sem1= std(sbeta2)./(sqrt(n_s)); % SEM beta-weights
% plot reference line going through zero on y-axis
plot([0 0],[0 n_x+1],'k-','LineWidth',2); hold on;
% plot error bars
o= 0; % offset
i_v=1:n_x; for i=1:length(i_v); plot([-1 +1],[i_v(i)-o i_v(i)-o],'k--','LineWidth',.1); hold on; end;
% me
i_v=1:7; for i=1:length(i_v); plot([mean1(i_v(i))-sem1(i_v(i)) mean1(i_v(i))+sem1(i_v(i))],[i_v(i)-o i_v(i)-o],'k-','LineWidth',2); hold on; end;
i_v=1:7; for i=1:length(i_v); plot(mean1(i_v(i)),i_v(i)-o,'ko','LineWidth',1,'MarkerFaceColor',me_col,'MarkerSize',12,'LineWidth',2); hold on; end;
% you
i_v=8:14; for i=1:length(i_v); plot([mean1(i_v(i))-sem1(i_v(i)) mean1(i_v(i))+sem1(i_v(i))],[i_v(i)-o i_v(i)-o],'k-','LineWidth',2); hold on; end;
i_v=8:14; for i=1:length(i_v); plot(mean1(i_v(i)),i_v(i)-o,'ko','LineWidth',1,'MarkerFaceColor',you_col,'MarkerSize',12,'LineWidth',2); hold on; end;
% general settings
set(gca,'YAxisLocation','left','FontSize',16,'LineWidth',2); % set font size for figure
box off; % only lines on x- and y-axis
% specify x-axis
set(gca,'XTick',[-.1:.1:.7]); % specify ticks (no need for labels as numeric)
xlabel('beta-weight','FontSize',36); % specify x-axis name
xlim([-.15 .65]); % set x-axis limits
% specify y-axis
set(gca,'YTick',1:n_x,'YTickLabel',predictor_names); % specify ticks and tick labels
ylabel('predictor','FontSize',36); % specify x-axis name
ylim([0 n_x+1]); % set y-axis limits
print('-djpeg','-r300',['Plots',filesep,'Figure6C']);
% statistics
[a2,b2,c2,d2]= ttest(sbeta2,0);
% print statistics
fprintf(['=================================================\n']);
fprintf(['HISTORY EFFECTS ON REVISED CONFIDENE\n']);
fprintf(['=================================================\n']);
fprintf(['p-vals: ' num2str(b2),'\n']);
fprintf(['t-vals: ' num2str(d2.tstat),'\n']);
fprintf(['=================================================\n']);

% VISUALISE HISTORY EFECTS ON REVISED CONFIDNECE
% calculate expected impact
% b1= +.2015;
% b2= +.0435;
% b3= -.0111;
% b4= +.0061;
% b5= -.0162;
% b6= +.0841;
% b7= +.0452;
muBeta= mean(sbeta2);
b1= muBeta(8);
b2= muBeta(9);
b3= muBeta(10);
b4= muBeta(11);
b5= muBeta(12);
b6= muBeta(13);
b7= muBeta(14);
c_v= [-3:1:3];
a_v= [-.5 .5];
for i_c1= 1:7
    for i_c2= 1:7
        for i_a= 1:2;
            impact(i_c1,i_c2,i_a)= b1*c_v(i_c1) + b2*a_v(i_a) + b3*c_v(i_c2) + b4*c_v(i_c2)*a_v(i_a) + b5*c_v(i_c1)*c_v(i_c2) + b6*c_v(i_c1)*a_v(i_a) + b7*c_v(i_c1)*c_v(i_c2)*a_v(i_a);
        end
    end
end
% figure
figure;
subplot(1,2,1);
imagesc(impact(:,:,1));
set(gca,'YDir','Normal');
set(gca,'FontSize',28,'LineWidth',2);
set(gca,'XTick',1:7,'XTickLabel',{'-3','-2','-1','0','1','2','3'});
set(gca,'YTick',1:7,'YTickLabel',{'-3','-2','-1','0','1','2','3'});
% ylabel('C_p_,_t','FontSize',36);
% xlabel('C_p_,_t_-_1','FontSize',36);
% title('A_p_,_t_-_1: incorrect','FontWeight','normal','FontSize',20);
% c=colorbar;
% set(c,'YTick',[],'LineWidth',2);
axis square;
% figure
subplot(1,2,2);
imagesc(impact(:,:,2));
set(gca,'YDir','Normal');
set(gca,'FontSize',28,'LineWidth',2);
set(gca,'XTick',1:7,'XTickLabel',{'-3','-2','-1','0','1','2','3'});
set(gca,'YTick',1:7,'YTickLabel',{'-3','-2','-1','0','1','2','3'});
% ylabel('C_p_,_t','FontSize',36);
% xlabel('C_p_,_t_-_1','FontSize',36);
% title('A_p_,_t_-_1: correct','FontWeight','normal','FontSize',20);
% c=colorbar;
% set(c,'YTick',[],'LineWidth',2);
axis square;
print('-djpeg','-r300',['Plots',filesep,'Figure6D']);