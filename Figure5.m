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
    
    % trial indices
    idx_advB= find(data.sbjID==i_sbj & data.advisor==1);
    idx_advG= find(data.sbjID==i_sbj & data.advisor==2);
    idx_both= find(data.sbjID==i_sbj);
    
    % performance: summary statistics
    % accuracy
    sacc1.advB(i_sbj) = mean(data.sbjAcc1(idx_advB));
    sacc1.advG(i_sbj) = mean(data.sbjAcc1(idx_advG));
    sacc1.both(i_sbj) = mean(data.sbjAcc1(idx_both));
    sacc2.advB(i_sbj) = mean(data.sbjAcc2(idx_advB));
    sacc2.advG(i_sbj) = mean(data.sbjAcc2(idx_advG));
    sacc2.both(i_sbj) = mean(data.sbjAcc2(idx_both));
    % confidence
    scon1.advB(i_sbj) = mean(abs(data.sbjCon1(idx_advB)));
    scon1.advG(i_sbj) = mean(abs(data.sbjCon1(idx_advG)));
    scon1.both(i_sbj) = mean(abs(data.sbjCon1(idx_both)));
    scon2.advB(i_sbj) = mean(abs(data.sbjCon2(idx_advG)));
    scon2.advG(i_sbj) = mean(abs(data.sbjCon2(idx_advG)));
    scon2.both(i_sbj) = mean(abs(data.sbjCon2(idx_both)));
    % score
    sscr1.advB(i_sbj) = mean(data.sbjScr1(idx_advB));
    sscr1.advG(i_sbj) = mean(data.sbjScr1(idx_advG));
    sscr1.both(i_sbj) = mean(data.sbjScr1(idx_both));
    sscr2.advB(i_sbj) = mean(data.sbjScr2(idx_advB));
    sscr2.advG(i_sbj) = mean(data.sbjScr2(idx_advG));
    sscr2.both(i_sbj) = mean(data.sbjScr2(idx_both));
    
    % social susceptibility: GLM
    % bad
    Y= zscore(data.sbjCon2(idx_advB))';
    X= [zscore(data.sbjCon1(idx_advB));
        zscore(data.sbjAcc1(idx_advB));
        zscore(data.advCon(idx_advB))]';
    B= glmfit(X,Y);
    sbeta.advB(i_sbj)= B(end);
    % good
    Y= zscore(data.sbjCon2(idx_advG))';
    X= [zscore(data.sbjCon1(idx_advG));
        zscore(data.sbjAcc1(idx_advG));
        zscore(data.advCon(idx_advG))]';
    B= glmfit(X,Y);
    sbeta.advG(i_sbj)= B(end);
    % both
    Y= zscore(data.sbjCon2(idx_both))';
    X= [zscore(data.sbjCon1(idx_both));
        zscore(data.sbjAcc1(idx_both));
        zscore(data.advCon(idx_both))]';
    B= glmfit(X,Y);
    sbeta.both(i_sbj)= B(end);
    
    % resolution
    % bad
    con1_v= abs(data.sbjCon1(idx_advB));
    acc1_v= data.sbjAcc1(idx_advB);
    sres1.advB(i_sbj)= mean(con1_v(acc1_v==1))-mean(con1_v(acc1_v==0));
    con2_v= abs(data.sbjCon2(idx_advB));
    acc2_v= data.sbjAcc2(idx_advB);
    sres2.advB(i_sbj)= mean(con2_v(acc1_v==1))-mean(con2_v(acc2_v==0));
    % bad
    con1_v= abs(data.sbjCon1(idx_advG));
    acc1_v= data.sbjAcc1(idx_advG);
    sres1.advG(i_sbj)= mean(con1_v(acc1_v==1))-mean(con1_v(acc1_v==0));
    con2_v= abs(data.sbjCon2(idx_advG));
    acc2_v= data.sbjAcc2(idx_advG);
    sres2.advG(i_sbj)= mean(con2_v(acc1_v==1))-mean(con2_v(acc2_v==0));
    % both
    con1_v= abs(data.sbjCon1(idx_both));
    acc1_v= data.sbjAcc1(idx_both);
    sres1.both(i_sbj)= mean(con1_v(acc1_v==1))-mean(con1_v(acc1_v==0));
    con2_v= abs(data.sbjCon2(idx_both));
    acc2_v= data.sbjAcc2(idx_both);
    sres2.both(i_sbj)= mean(con2_v(acc1_v==1))-mean(con2_v(acc2_v==0));
    
    % calibration
    % bad
    con1_v= abs(data.sbjCon1(idx_advB));
    acc1_v= data.sbjAcc1(idx_advB);
    con1_u= unique(con1_v);
    tmp= []; 
    for i= 1:length(con1_u); 
        tmp(i)= mean(con1_v(con1_v==con1_u(i))) - mean(acc1_v(con1_v==con1_u(i)));
    end
    scal1.advB(i_sbj)= sum(tmp);
    con2_v= abs(data.sbjCon2(idx_advB));
    acc2_v= data.sbjAcc2(idx_advB);
    con2_u= unique(con2_v);
    tmp= []; 
    for i= 1:length(con2_u); 
        tmp(i)= mean(con2_v(con2_v==con2_u(i))) - mean(acc2_v(con2_v==con2_u(i)));
    end
    scal2.advB(i_sbj)= sum(tmp);
    % bad
    con1_v= abs(data.sbjCon1(idx_advG));
    acc1_v= data.sbjAcc1(idx_advG);
    con1_u= unique(con1_v);
    tmp= []; 
    for i= 1:length(con1_u); 
        tmp(i)= mean(con1_v(con1_v==con1_u(i))) - mean(acc1_v(con1_v==con1_u(i)));
    end
    scal1.advG(i_sbj)= sum(tmp);
    con2_v= abs(data.sbjCon2(idx_advG));
    acc2_v= data.sbjAcc2(idx_advG);
    con2_u= unique(con2_v);
    tmp= []; 
    for i= 1:length(con2_u); 
        tmp(i)= mean(con2_v(con2_v==con2_u(i))) - mean(acc2_v(con2_v==con2_u(i)));
    end
    scal2.advG(i_sbj)= sum(tmp);
    % bad
    con1_v= abs(data.sbjCon1(idx_both));
    acc1_v= data.sbjAcc1(idx_both);
    con1_u= unique(con1_v);
    tmp= []; 
    for i= 1:length(con1_u); 
        tmp(i)= mean(con1_v(con1_v==con1_u(i))) - mean(acc1_v(con1_v==con1_u(i)));
    end
    scal1.both(i_sbj)= sum(tmp);
    con2_v= abs(data.sbjCon2(idx_both));
    acc2_v= data.sbjAcc2(idx_both);
    con2_u= unique(con2_v);
    tmp= []; 
    for i= 1:length(con2_u); 
        tmp(i)= mean(con2_v(con2_v==con2_u(i))) - mean(acc2_v(con2_v==con2_u(i)));
    end
    scal2.both(i_sbj)= sum(tmp);
    
    % confidence by accuracy for group-level plot
    cxa.adv{1}.rsp{1}.acc{1}(i_sbj) = mean(abs(data.sbjCon1(data.sbjAcc1==0 & data.advisor==1)));
    cxa.adv{1}.rsp{2}.acc{1}(i_sbj) = mean(abs(data.sbjCon2(data.sbjAcc2==0 & data.advisor==1)));
    cxa.adv{1}.rsp{1}.acc{2}(i_sbj) = mean(abs(data.sbjCon1(data.sbjAcc2==1 & data.advisor==1)));
    cxa.adv{1}.rsp{2}.acc{2}(i_sbj) = mean(abs(data.sbjCon2(data.sbjAcc2==1 & data.advisor==1)));
    cxa.adv{2}.rsp{1}.acc{1}(i_sbj) = mean(abs(data.sbjCon1(data.sbjAcc1==0 & data.advisor==2)));
    cxa.adv{2}.rsp{2}.acc{1}(i_sbj) = mean(abs(data.sbjCon2(data.sbjAcc2==0 & data.advisor==2)));
    cxa.adv{2}.rsp{1}.acc{2}(i_sbj) = mean(abs(data.sbjCon1(data.sbjAcc2==1 & data.advisor==2)));
    cxa.adv{2}.rsp{2}.acc{2}(i_sbj) = mean(abs(data.sbjCon2(data.sbjAcc2==1 & data.advisor==2)));
    
    
end

%% -----------------------------------------------------------------------
%% FIGURE
%% -----------------------------------------------------------------------

% SUSCEPTIILITY ~ CALIBRATION
% figure
figz=figure('color',[1 1 1]);
% colours
adv_col= [1.00 0.85 0.34; ... % bad B: yellow
          0.10 0.60 1.00];    % good G: blue
% individual statistics
Y= sbeta.both';
X= [sres1.both; scal1.both; sacc1.both]';
[B STATS]= robustfit(zscore(X),zscore(Y));
% plot individual statistics
plot([-3 3],[0 0],'k--'); hold on;
plot([0 0],[-3 3],'k--'); hold on;
hold on;
plot(zscore(scal1.both),zscore(Y),'ko','MarkerFaceColor',[255,140,0]./255,'MarkerSize',12);
plot(zscore(scal1.both),B(1)+B(3)*zscore(scal1.both),'k-','LineWidth',2); hold on; 
% specify axes
xlim([-3 3]);
ylim([-3 3]);
set(gca,'XTick',-3:1:3);
set(gca,'YTick',-3:1:3);
set(gca,'LineWidth',2,'FontSize',28);
box off;
xlabel('calibration [z-scored]','FontSize',36);
ylabel('susceptibility [z-scored]','FontSize',36);
% print figure
print('-djpeg','-r300',['Plots',filesep,'Figure6A']);

% SUSCEPTIILITY ~ RESOLUTION
% figure
figz=figure('color',[1 1 1]);
% colours
adv_col= [1.00 0.85 0.34; ... % bad B: yellow
          0.10 0.60 1.00];    % good G: blue
% individual statistics
Y= sbeta.both';
X= [sres1.both; scal1.both; sacc1.both]';
[B STATS]= robustfit(zscore(X),zscore(Y));
% plot individual statistics
plot([-3 3],[0 0],'k--'); hold on;
plot([0 0],[-3 3],'k--'); hold on;
hold on;
plot(zscore(sres1.both),zscore(Y),'ko','MarkerFaceColor',[255,140,0]./255,'MarkerSize',12);
plot(zscore(sres1.both),B(1)+B(2)*zscore(sres1.both),'k-','LineWidth',2); hold on; 
% specify axes
xlim([-3 3]);
ylim([-3 3]);
set(gca,'XTick',-3:1:3);
set(gca,'YTick',-3:1:3);
set(gca,'LineWidth',2,'FontSize',28);
box off;
xlabel('resolution [z-scored]','FontSize',36);
ylabel('susceptibility [z-scored]','FontSize',36);
% print figure
print('-djpeg','-r300',['Plots',filesep,'Figure6B']);

% SUSCEPTIILITY ~ CALIBRATION
% figure
figz=figure('color',[1 1 1]);
% colours
adv_col= [1.00 0.85 0.34; ... % bad B: yellow
          0.10 0.60 1.00];    % good G: blue
% individual statistics
Y= [sscr2.both-sscr1.both]';
X= [sres1.both; scal1.both; sacc1.both]';
[B STATS]= robustfit(zscore(X),zscore(Y));
% plot individual statistics
plot([-3 3],[0 0],'k--'); hold on;
plot([0 0],[-3 3],'k--'); hold on;
hold on;
plot(zscore(scal1.both),zscore(Y),'ko','MarkerFaceColor',[255,140,0]./255,'MarkerSize',12);
plot(zscore(scal1.both),B(1)+B(3)*zscore(scal1.both),'k-','LineWidth',2); hold on; 
% specify axes
xlim([-3 3]);
ylim([-3 3]);
set(gca,'XTick',-3:1:3);
set(gca,'YTick',-3:1:3);
set(gca,'LineWidth',2,'FontSize',28);
box off;
xlabel('calibration [z-scored]','FontSize',36);
ylabel('benefit [z-scored]','FontSize',36);
% print figure
print('-djpeg','-r300',['Plots',filesep,'Figure6C']);

% SUSCEPTIILITY ~ RESOLUTION
% figure
figz=figure('color',[1 1 1]);
% colours
adv_col= [1.00 0.85 0.34; ... % bad B: yellow
          0.10 0.60 1.00];    % good G: blue
% individual statistics
Y= [sscr2.both-sscr1.both]';
X= [sres1.both; scal1.both; sacc1.both]';
[B STATS]= robustfit(zscore(X),zscore(Y));
% plot individual statistics
plot([-3 3],[0 0],'k--'); hold on;
plot([0 0],[-3 3],'k--'); hold on;
hold on;
plot(zscore(sres1.both),zscore(Y),'ko','MarkerFaceColor',[255,140,0]./255,'MarkerSize',12);
plot(zscore(sres1.both),B(1)+B(2)*zscore(sres1.both),'k-','LineWidth',2); hold on; 
% specify axes
xlim([-3 3]);
ylim([-3 3]);
set(gca,'XTick',-3:1:3);
set(gca,'YTick',-3:1:3);
set(gca,'LineWidth',2,'FontSize',28);
box off;
xlabel('resolution [z-scored]','FontSize',36);
ylabel('benefit [z-scored]','FontSize',36);
% print figure
print('-djpeg','-r300',['Plots',filesep,'Figure6D']);