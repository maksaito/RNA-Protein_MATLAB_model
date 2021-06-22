% RNA_protein_model_062221.m 
% Version  06/22/2021
% Model by Mak Saito, WHOI 2018-2021 for MATLAB
% Simulating RNA and protein production, inventories, decay, and influence on fold-change
% Used in Walworth, Saito et al., 2021 submitted

clear

% CONSTANTS ---------------------------------------------
% time steps: 3 days in minute interals (1440 minutes per day *3)
time=[1:4320]'; %'

% *** SIGNAL parameters *** 0-1 scale
	% s0 = basal expression
	s0=0.01; % set at 1 percent of cells (non-zero) to allow ratios
	% s1 = expression at event 1
	s1=0.2; % 20 percent of cells respond to stimuli

	% alternate values for basal expression simulations
	s0alt=0.2;
	%s1alt=0.3;		
	s1alt=0.4;		

	% Duration of each signal level
	
	% Select ramp by removing comment symbol %
	% Short ramp:
	tstart=600; tplateau1=630; tplateau2=720; tend=750;
	
	% Slow ramp:
	%tstart=600; tplateau1=700; tplateau2=800; tend=1200;
	
	% Long slow ramp:
	tstart=600; tplateau1=650; tplateau2=700; tend=1200;

	% Long slow ramp:
	%tstart=600; tplateau1=650; tplateau2=700; tend=2400;

	% Alternate 1 minute signal for testing:
	%tstart=630;
	%tplateau1=630;
	%tplateau2=631;
	%tend=631;

	% calculate linear change in signal between tstart & tplateau1: amount added per time step
	sadd = ((s1-s0)/((tplateau1-tstart)));
	saddalt = ((s1alt-s0alt)/((tplateau1-tstart)));
	% calculate linear change in signal between tplateau2 & tend: amount added per time step
	sminus = ((s1-s0)/((tplateau2-tend)));
	sminusalt = ((s1alt-s0alt)/((tplateau2-tend)));

% *** Production of RNA and Protein *** 
	krna = 100; % production of a RNA transcript per 100 cells per time step
	kprotein = 3.5; % tuned using 1 minute signal to be ~6-7 proteins per transcript production based on Moran et al. 2013

% *** Degradation of RNA and Protein ***
	% krnaloss = 0.462; % (min-1) = 1.5 min half life Moran et al 2013 (from bionumbers)
	% kprotloss=0.000963 %9.63e-4; % (min-1) = 12h half life Moran et al., 2013 
	% half life to rate constant conversion: A = Ao*e^(kt); A = 0.5Ao; ln(0.5)= -0.69; t1/2 = -0.69 / k; careful of sign error
	krnaloss = 0.102; % (min-1) 6.8 min half life Selinger et al., 2003 Genome Res. E coli
	kprotloss=0.0017; %  6.9 - 34 h range half lives (k = 0.1 - 0.022 h-1; 0.0017- 3.7e-4 m-1) Pratt et al., 2002 Mol Cell Proteomics

	
% *** Delay between signal and RNA ***
	drna=1; % delay constants in minutes
	% delay between transcription and translation
	dprotein=180; % delay constants in minutes	

% initialize arrays
	% signal to s0
	signal(:,1)=ones(4320,1).*s0;
	signal(:,2)=ones(4320,1).*s0alt;
	binsignal=zeros(4320,2);
	RNA=zeros(4320,2);
	protein=zeros(4320,2);
	protein(1:20,2)=1800; %initial protein concentration (first 20min due to time steps)
	binRNA=zeros(4320,2);
	binprotein=zeros(4320,2);

% select time zero point for fold-change calculations
	% one hour before signal initiates
	t0=9*60; % 9 hours
%---------------------------------------------------------------------------------------

% *** Calculate Signal ***
% First panels: Initial signal = 0 
% Middle Panels: Initial signal > zero (s0alt), increases to s1alt 
for i = 1:4320
	if time(i) > tstart
		signal(i,1) = signal(i-1,1)+sadd;
		signal(i,2) = signal(i-1,2)+saddalt;
	if time(i) > tplateau1
		signal(i,1)=s1;
		signal(i,2)=s1alt;
	if time(i) > tplateau2
		signal(i,1) = signal(i-1,1) + sminus;
		signal(i,2) = signal(i-1,2) + sminusalt;
	if time(i) > tend
		signal(i,1) = s0;
		signal(i,2) = s0alt;
	if time(i) > tstart+1440
		signal(i,1) = signal(i-1,1)+sadd;
		signal(i,2) = signal(i-1,2)+saddalt;
	if time(i) > tplateau1+1440
		signal(i,1)=s1;
		signal(i,2)=s1alt;
	if time(i) > tplateau2+1440
		signal(i,1) = signal(i-1,1) + sminus;
		signal(i,2) = signal(i-1,2) + sminusalt;
	if time(i) > tend+1440
		signal(i,1) = s0;
		signal(i,2) = s0alt;
	if time(i) > tstart+2880
		signal(i,1) = signal(i-1,1)+sadd;
		signal(i,2) = signal(i-1,2)+saddalt;
	if time(i) > tplateau1+2880
		signal(i,1)=s1;
		signal(i,2)=s1alt;
	if time(i) > tplateau2+2880
		signal(i,1) = signal(i-1,1) + sminus;
		signal(i,2) = signal(i-1,2) + sminusalt;
	if time(i) > tend+2880
		signal(i,1) = s0;
		signal(i,2) = s0alt;
	end
	end
	end
	end
	end
	end
	end
	end
	end
	end
	end
	end
end

% *** RNA Inventory Calculations ***
% calculate RNA (1) no basal expression (2) basal expression, using offset drna 
for i = (drna+1):4320  % start index at i >=2
	RNA(i,1) = RNA(i-1,1) + krna*signal(i-drna,1)-krnaloss*RNA(i-1,1);
	RNA(i,2) = RNA(i-1,2) + krna*signal(i-drna,2)-krnaloss*RNA(i-1,2);
end

% *** Protein Inventory Calculations ***
% calculate Protein (1) no basal expression (2) basal expression, using offset dprotein 
% prior protein time step, plus production from RNA (with time step 20min), minus protein loss prior step
for i = (dprotein+1):4320
	protein(i,1) = protein(i-1,1) + kprotein*RNA(i-dprotein,1) - kprotloss*protein(i-1,1);
	protein(i,2) = protein(i-1,2) + kprotein*RNA(i-dprotein,2) - kprotloss*protein(i-1,2);
end

% Calculate fold change relative to t0
	RNAfold_s0_zero = RNA(:,1)./RNA(t0,1);
	Proteinfold_s0_zero = protein(:,1)./protein(t0,1);
	RNAfold_s0_notzero = RNA(:,2)./RNA(t0,2);
	Proteinfold_s0_notzero = protein(:,2)./protein(t0,2);

% Calculate fold change each day 
% start by settting new arrays divided by t0
	RNAfold_s0_zero_byday = RNA(:,1)./RNA(t0,1);
	Proteinfold_s0_zero_byday = protein(:,1)./protein(t0,1);
	RNAfold_s0_notzero_byday = RNA(:,2)./RNA(t0,2);
	Proteinfold_s0_notzero_byday = protein(:,2)./protein(t0,2);
% rewrite day 2
for i = 1441:2880
	RNAfold_s0_zero_byday(i,1) = RNA(i,1)/RNA(t0+1440,1);
	Proteinfold_s0_zero_byday(i,1) = protein(i,1)/protein(t0+1440,1);
	RNAfold_s0_notzero_byday(i,1) = RNA(i,2)/RNA(t0+1440,2);
	Proteinfold_s0_notzero_byday(i,1) = protein(i,2)/protein(t0+1440,2);
end
% rewrite day 3
for i = 2881:4320
	RNAfold_s0_zero_byday(i,1) = RNA(i,1)/RNA(t0+2880,1);
	Proteinfold_s0_zero_byday(i,1) = protein(i,1)/protein(t0+2880,1);
	RNAfold_s0_notzero_byday(i,1) = RNA(i,2)/RNA(t0+1440,2);
	Proteinfold_s0_notzero_byday(i,1) = protein(i,2)/protein(t0+2880,2);
end

% --- FIGURES -------------------------------------------------------------

% choose time window (i= end) to plot up to 4320 minutes (3 days)
%i=1440
i=4320;
figure(1)
clf
subplot(321)
plot(time(1:i), signal(1:i,1),'k-')
%axis([0 20 0 1.2])
title('1% Basal Signal')
ylabel('Signal')
subplot(323)
axis([0 20 0 2])
plot(time(1:i),RNA(1:i),'r-')
ylabel('Transcripts/100 cells')
subplot(325)
plot(time(1:i),protein(1:i),'b-')
%axis([0 20 0 70])
ylabel('Proteins/100 cells')
xlabel('Time (minutes)')
subplot(322)
plot(time(1:i), signal(1:i,2),'k-')
%axis([0 20 0 1.2])
title('20% Basal Signal')
ylabel('Signal')
subplot(324)
plot(time(1:i),RNA(1:i,2),'r-')
%axis([0 20 0 2])
ylabel('Transcripts/100 cells')
subplot(326)
plot(time(1:i),protein(1:i,2),'b-')
%axis([0 20 0 70])
ylabel('Proteins/100 cells')
xlabel('Time (minutes)')
print('rnap_fig1', '-dpng', '-r600'); %


figure(2)
clf
subplot(221)
plot(time(1:i), RNAfold_s0_zero(1:i), 'r-')
xlabel('minutes')
ylabel('RNA fold change')
title('1% Basal Signal (cells/min)')
subplot(223)
plot(time(1:i), Proteinfold_s0_zero(1:i), 'b-')
xlabel('minutes')
ylabel('Protein fold change')
subplot(222)
plot(time(1:i), RNAfold_s0_notzero(1:i), 'r-'); hold on
plot(time(1:i), RNAfold_s0_notzero_byday(1:i), 'k.-')
xlabel('minutes')
ylabel('RNA fold change')
title('20% Basal Signal (cells/min)')
subplot(224)
plot(time(1:i), Proteinfold_s0_notzero(1:i), 'b-'); hold on
plot(time(1:i), Proteinfold_s0_notzero_byday(1:i), 'k-')
xlabel('minutes')
ylabel('Protein fold change')
legend('/ day 0','/ each day','location','northwest') 
print('rnap_fig2', '-dpng', '-r600'); %

figure(3)
clf
xline=[-10,0; 50, 0]; 
yline=[0, -10; 0, 15];
onetoone=[-10; 50];
subplot(211)
plot(RNAfold_s0_zero,Proteinfold_s0_zero,'ob'); hold on
plot(xline(:,1), xline(:,2),'k-');
plot(yline(:,1), yline(:,2),'k-');
plot(onetoone(:,1), onetoone(:,1),'-.')
axis([-10 50 -10 50])
xlabel('RNA fold change')
ylabel('Protein fold change')
title('No pre-existing inventory')
subplot(212)
plot(RNAfold_s0_notzero,Proteinfold_s0_notzero,'ob'); hold on
plot(xline(:,1), xline(:,2),'k-');
plot(yline(:,1), yline(:,2),'k-');
plot(onetoone(:,1), onetoone(:,1),'-.')
axis([-10 50 -10 15])
xlabel('RNA fold change')
ylabel('Protein fold change')
title('Pre-existing inventory')

figure(4)
clf
xline=[-10,0; 50, 0]; 
yline=[0, -10; 0, 15];
onetoone=[-2; 7];

subplot(231)
plot(log(RNAfold_s0_zero(1:1440)),log(Proteinfold_s0_zero(1:1440)),'b.-'); hold on
plot(log(RNAfold_s0_zero(1441:2880)),log(Proteinfold_s0_zero(1441:2880)),'r.-'); hold on
plot(log(RNAfold_s0_zero(2881:4320)),log(Proteinfold_s0_zero(2881:4320)),'g.-'); hold on
plot(log(RNAfold_s0_zero(963)),log(Proteinfold_s0_zero(963)),'m*')
plot(log(RNAfold_s0_zero(2060)),log(Proteinfold_s0_zero(2060)),'kd')
plot(xline(:,1), xline(:,2),'k-');
plot(yline(:,1), yline(:,2),'k-');
plot(onetoone(:,1), onetoone(:,1),'-.')
axis([-2 10 -1 5])
xlabel('Log2 RNA fold change')
ylabel('Log2 Protein fold change')
%title('No pre-existing inventory; relative to day 1')
title('1% basal expression - relative to day 1')
legend('day 1','day 2', 'day 3','min 963', 'min 2060')

subplot(232)
%plot(log(RNAfold_s0_notzero),log(Proteinfold_s0_notzero),'r-+'); hold on
plot(log(RNAfold_s0_notzero_byday(1:1440)),log(Proteinfold_s0_notzero_byday(1:1440)),'b.-'); hold on
plot(log(RNAfold_s0_notzero_byday(1441:2880)),log(Proteinfold_s0_notzero_byday(1441:2880)),'r.-'); hold on
plot(log(RNAfold_s0_notzero_byday(2881:4320)),log(Proteinfold_s0_notzero_byday(2881:4320)),'g.-'); hold on
plot(log(RNAfold_s0_notzero_byday(860)),log(Proteinfold_s0_notzero_byday(860)),'m*')
plot(log(RNAfold_s0_notzero_byday(2500)),log(Proteinfold_s0_notzero_byday(2500)),'kd')
plot(xline(:,1), xline(:,2),'k-');
plot(yline(:,1), yline(:,2),'k-');
plot(onetoone(:,1), onetoone(:,1),'-.')
axis([-2 10 -1 5])
%axis([-2 2 -1 1])
xlabel('Log2 RNA fold change')
ylabel('Log2 Protein fold change')
title('20% Basal expression')
legend('day 1','day 2', 'day 3','min 860', 'min 2500')

subplot(234)
plot(log(RNAfold_s0_zero_byday(1:1440)),log(Proteinfold_s0_zero_byday(1:1440)),'b.-'); hold on
plot(log(RNAfold_s0_zero_byday(1441:2880)),log(Proteinfold_s0_zero_byday(1441:2880)),'r.-'); hold on
plot(log(RNAfold_s0_zero_byday(2881:4320)),log(Proteinfold_s0_zero_byday(2881:4320)),'g.-'); hold on
plot(log(RNAfold_s0_zero_byday(963)),log(Proteinfold_s0_zero_byday(963)),'m*')
plot(log(RNAfold_s0_zero_byday(2600)),log(Proteinfold_s0_zero_byday(2600)),'kd')
plot(xline(:,1), xline(:,2),'k-');
plot(yline(:,1), yline(:,2),'k-');
plot(onetoone(:,1), onetoone(:,1),'-.')
axis([-2 10 -1 5])
xlabel('Log2 RNA fold change')
ylabel('Log2 Protein fold change')
%title('No pre-existing inventory ; relative to each day')
title('1% basal expression - relative to each day')
legend('day 1','day 2', 'day 3','min 963', 'min 2600')

subplot(235)
%plot(log(RNAfold_s0_notzero(1:1440)),log(Proteinfold_s0_notzero(1:1440)),'k-+'); hold on
plot(log(RNAfold_s0_notzero_byday(1:1440)),log(Proteinfold_s0_notzero_byday(1:1440)),'b.-'); hold on
plot(log(RNAfold_s0_notzero_byday(1441:2880)),log(Proteinfold_s0_notzero_byday(1441:2880)),'r.-'); hold on
plot(log(RNAfold_s0_notzero_byday(2881:4320)),log(Proteinfold_s0_notzero_byday(2881:4320)),'g.-'); hold on
plot(log(RNAfold_s0_notzero_byday(860)),log(Proteinfold_s0_notzero_byday(860)),'m*')
plot(log(RNAfold_s0_notzero_byday(2500)),log(Proteinfold_s0_notzero_byday(2500)),'kd')
plot(xline(:,1), xline(:,2),'k-');
plot(yline(:,1), yline(:,2),'k-');
plot(onetoone(:,1), onetoone(:,1),'-.')
% cross 1:1 line
axis([-2 2 -1 1])
xlabel('Log2 RNA fold change')
ylabel('Log2 Protein fold change')
%title('Pre-existing inventory; zoomed in')
title('20% Basal expression - expanded axes')
legend('day 1','day 2', 'day 3', 'min 860', 'min 2500','Location','Northwest')

subplot(233)
foldmatrix=[(RNAfold_s0_zero(963)),(Proteinfold_s0_zero(963)); ...
	 (RNAfold_s0_zero(2060)),(Proteinfold_s0_zero(2060)); ...
	 (RNAfold_s0_notzero_byday(860)),(Proteinfold_s0_notzero_byday(860))];
bar(foldmatrix,0.75,'grouped'); hold on
xline=[0 4]; yline=[2 2];
plot(xline, yline,'k:')
ylabel('Fold change')
str = {'1% day1';'1% day2'; '20%'};
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
legend('RNA','Protein','Two-fold')
title('Fold change relative to day 1')

subplot(236)
foldmatrix_byday=[(RNAfold_s0_zero(963)),(Proteinfold_s0_zero(963)); ...
	(RNAfold_s0_zero_byday(2600)),(Proteinfold_s0_zero_byday(2600)); ...
	(RNAfold_s0_notzero_byday(860)),(Proteinfold_s0_notzero_byday(860))];
bar(foldmatrix_byday,0.75,'grouped'); hold on
plot(xline, yline,'k:')
ylabel('Fold change')
str = {'1% day1';'1% day2'; '20%'};
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
legend('RNA','Protein','Two-fold')
title('Fold change relative to each day')
print('rnap_fig4', '-dpng', '-r600'); %


