% Matlab code for:

% Alexander Heyde, David Rohde, Cameron S. McAlpine, Shuang Zhang, Friedrich F. Hoyer, Jeffrey M. Gerold, David Cheek, Yoshiko Iwamoto, Maximilian J. Schloss, Katrien Vandoorne, Oriol Iborra - Egea, Christian Muñoz - Guijosa, Antoni Bayes - Genis, Johannes G. Reiter, Morgan Craig, Filip K. Swirski, Matthias Nahrendorf, Martin A. Nowak, Kamila Naxerova. 
% Increased stem cell proliferation in atherosclerosis accelerates clonal hematopoiesis.

%% Figure 3F

% Table S2
% All time measurements in years
b = 365/28;   % Average baseline HSC proliferation rate
T = 40;       % Age of onset for elevated HSC proliferation rate
F = .02;      % Detection frequency (minimum VAF)
ts = 50;      % Mean age of baseline VAF data
fs = .0051;   % Frequency of the largest driver clone at age ts
s = .003;     % Selective effect of the largest driver clone
N = 10^5;     % Number of HSCs

t1 = 70;      % Age at sequencing
M = 15000;    % Number of realizations
dt = .01;     % Time step for simuation

% Set distribution of proliferation rates for control cohort
r = [0.4883	1.8066	0.9861	0.9861	0.9151	1.0966	0.2832	0.6453	1.9959	0.7968]; %from Fig 2C
R = normrnd(mean(r), std(r), 1, M);
x = 2.^(randn(1,M)-7.47033)+zeros(1,M); % Initialization

% Simulate driver clone expansion in control cohort
for t = dt:dt:(b*(t1-T))
    x = x + (s*x.*(1-x) + sqrt(x.*(1-x)/N).*randn(1,M)).*R*dt;
    x = min(1,max(0,x));
end
Px = sum(x/2 > F)/M;
X = mean(x/2);

% Plot driver clone VAF in control cohort
figure('Position',[500 350 650 350])
edges = -4:.1:0;
histogram(x/2,10.^edges,'Normalization','probability')
set(gca,'xscale','log')
xlim(10.^[-3 0])
xlabel('Driver VAF')
ylabel('Fraction of patients')
hold on

% Set distribution of proliferation rates for disease cohort
r = [1.9328	2.2878	1.6961	3.8025	1.4595	3.3843	1.7908	3.7393	1.9170	2.8637]; %from Fig 2C
Ry = normrnd(mean(r), std(r), 1, M);
y = 2.^(randn(1,M)-7.47033)+zeros(1,M); % Initialization

% Simulate driver clone expansion in disease cohort
for t = dt:dt:(b*(t1-T))
    y = y + (s*y.*(1-y) + sqrt(y.*(1-y)/N).*randn(1,M)).*Ry*dt;
    y = min(1,max(0,y));
end
Py = sum(y/2 > F)/M;
Y = mean(y/2);

% Plot driver clone VAF in disease cohort
histogram(y/2,10.^edges,'Normalization','probability')
set(gca,'xscale','log')
xlim([.0001 1])
ylim([0 .12])
hold on
line([F F],[0 1],'Color','red','LineStyle','--')
