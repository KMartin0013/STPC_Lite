% make the coefficients for simulation
% note that this is a case for ocean usage

% clc;clear;

%% In the Boundary
%1
Coeffin(1).name='in';
Coeffin(1).Coef_trend=1/365.25;
Coeffin(1).Coef_A=2;
Coeffin(1).Coef_w=2*pi/365.25;
Coeffin(1).Coef_fi=pi/6;
Coeffin(1).Coef_Eps=1;

%2 similar to 1 with smaller noise
Coeffin(2).name='in';
Coeffin(2).Coef_trend=1/365.25;
Coeffin(2).Coef_A=2;
Coeffin(2).Coef_w=2*pi/365.25;
Coeffin(2).Coef_fi=pi/6;
Coeffin(2).Coef_Eps=0.1;

%3 only seasonal without noise
Coeffin(3).name='in';
Coeffin(3).Coef_trend=0;
Coeffin(3).Coef_A=0;
Coeffin(3).Coef_w=2*pi/365.25;
Coeffin(3).Coef_fi=pi/6;
Coeffin(3).Coef_Eps=0;

%4 trend and seasonal with noise
Coeffin(4).name='in';
Coeffin(4).Coef_trend=2/365.25;
Coeffin(4).Coef_A=2;
Coeffin(4).Coef_w=2*pi/365.25;
Coeffin(4).Coef_fi=pi/6;
Coeffin(4).Coef_Eps=1;

%5 trend (higher) and seasonal with noise
Coeffin(5).name='in';
Coeffin(5).Coef_trend=100/365.25;
Coeffin(5).Coef_A=2;
Coeffin(5).Coef_w=2*pi/365.25;
Coeffin(5).Coef_fi=pi/6;
Coeffin(5).Coef_Eps=1;

%6 only trend with noise
Coeffin(6).name='in';
Coeffin(6).Coef_trend=1/365.25;
Coeffin(6).Coef_A=0;
Coeffin(6).Coef_w=0;
Coeffin(6).Coef_fi=0;
Coeffin(6).Coef_Eps=1;

%7 similar to 2 but with decreasing trend
Coeffin(7).name='in';
Coeffin(7).Coef_trend=-1/365.25;
Coeffin(7).Coef_A=2;
Coeffin(7).Coef_w=2*pi/365.25;
Coeffin(7).Coef_fi=pi/6;
Coeffin(7).Coef_Eps=0.1;

%8 similar to 2 but with larger noise
Coeffin(8).name='in';
Coeffin(8).Coef_trend=1/365.25;
Coeffin(8).Coef_A=2;
Coeffin(8).Coef_w=2*pi/365.25;
Coeffin(8).Coef_fi=pi/6;
Coeffin(8).Coef_Eps=1;

%% Out of the Boundary
% 1 no signal of environment
Coeffout(1).name='out';
Coeffout(1).Coef_trend=0;
Coeffout(1).Coef_A=0;
Coeffout(1).Coef_w=0;
Coeffout(1).Coef_fi=0;
Coeffout(1).Coef_Eps=0.1;

% only seasonal with noise
Coeffout(2).name='out';
Coeffout(2).Coef_trend=0;
Coeffout(2).Coef_A=0;
Coeffout(2).Coef_w=2*pi/365.25;
Coeffout(2).Coef_fi=pi/6;
Coeffout(2).Coef_Eps=1;

% only seasonal without noise
Coeffout(3).name='out';
Coeffout(3).Coef_trend=0;
Coeffout(3).Coef_A=0;
Coeffout(3).Coef_w=2*pi/365.25;
Coeffout(3).Coef_fi=pi/6;
Coeffout(3).Coef_Eps=0;

% only trend with noise
Coeffout(4).name='out';
Coeffout(4).Coef_trend=1/365.25;
Coeffout(4).Coef_A=0;
Coeffout(4).Coef_w=0;
Coeffout(4).Coef_fi=0;
Coeffout(4).Coef_Eps=1;

% (5:10) corresponding to In 1 with changing Trend of environment
% Out 5 same to In 1
for i=1:5
Coeffout(i+4).name='out';
Coeffout(i+4).Coef_trend=1/365.25*i;
Coeffout(i+4).Coef_A=2;
Coeffout(i+4).Coef_w=2*pi/365.25;
Coeffout(i+4).Coef_fi=pi/6;
Coeffout(i+4).Coef_Eps=1;
end

Coeffout(10).name='out';
Coeffout(10).Coef_trend=1/365.25*0.5;
Coeffout(10).Coef_A=2;
Coeffout(10).Coef_w=2*pi/365.25;
Coeffout(10).Coef_fi=pi/6;
Coeffout(10).Coef_Eps=1;

% (11:16) corresponding to In 1 with changing Seasonal amplitude of environment
% Out 11 same to In 1
for i=1:5
Coeffout(i+10).name='out';
Coeffout(i+10).Coef_trend=1/365.25;
Coeffout(i+10).Coef_A=2*i;
Coeffout(i+10).Coef_w=2*pi/365.25;
Coeffout(i+10).Coef_fi=pi/6;
Coeffout(i+10).Coef_Eps=1;
end

Coeffout(16).name='out';
Coeffout(16).Coef_trend=1/365.25;
Coeffout(16).Coef_A=1;
Coeffout(16).Coef_w=2*pi/365.25;
Coeffout(16).Coef_fi=pi/6;
Coeffout(16).Coef_Eps=1;

% (17:21) corresponding to In 1 with changing Seasonal phase of environment
% Out 17 same to In 1
for i=1:3
Coeffout(i+16).name='out';
Coeffout(i+16).Coef_trend=1/365.25;
Coeffout(i+16).Coef_A=2;
Coeffout(i+16).Coef_w=2*pi/365.25;
Coeffout(i+16).Coef_fi=pi/6+pi/6*(i-1);
Coeffout(i+16).Coef_Eps=1;
end

for i=1:2
Coeffout(i+19).name='out';
Coeffout(i+19).Coef_trend=1/365.25;
Coeffout(i+19).Coef_A=2;
Coeffout(i+19).Coef_w=2*pi/365.25;
Coeffout(i+19).Coef_fi=pi/6-pi/6*i;
Coeffout(i+19).Coef_Eps=1;
end

Coeffout(22).name='out';
Coeffout(22).Coef_trend=1/365.25*0.25;
Coeffout(22).Coef_A=2;
Coeffout(22).Coef_w=2*pi/365.25;
Coeffout(22).Coef_fi=pi/6;
Coeffout(22).Coef_Eps=1;

Coeffout(23).name='out';
Coeffout(23).Coef_trend=1/365.25;
Coeffout(23).Coef_A=0.5;
Coeffout(23).Coef_w=2*pi/365.25;
Coeffout(23).Coef_fi=pi/6;
Coeffout(23).Coef_Eps=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (24:29) corresponding to In 2 with changing Trend of environment (for sensitivity analysis)
% Out 5 same to In 1
for i=1:5
Coeffout(i+23).name='out';
Coeffout(i+23).Coef_trend=1/365.25*i;
Coeffout(i+23).Coef_A=2;
Coeffout(i+23).Coef_w=2*pi/365.25;
Coeffout(i+23).Coef_fi=pi/6;
Coeffout(i+23).Coef_Eps=0.1;
end

Coeffout(29).name='out';
Coeffout(29).Coef_trend=1/365.25*0.5;
Coeffout(29).Coef_A=2;
Coeffout(29).Coef_w=2*pi/365.25;
Coeffout(29).Coef_fi=pi/6;
Coeffout(29).Coef_Eps=0.1;

% (30:35) corresponding to In 2 with changing Seasonal amplitude of
% environment (for sensitivity analysis)
% Out 11 same to In 1
for i=1:5
Coeffout(i+29).name='out';
Coeffout(i+29).Coef_trend=1/365.25;
Coeffout(i+29).Coef_A=2*i;
Coeffout(i+29).Coef_w=2*pi/365.25;
Coeffout(i+29).Coef_fi=pi/6;
Coeffout(i+29).Coef_Eps=0.1;
end

Coeffout(35).name='out';
Coeffout(35).Coef_trend=1/365.25;
Coeffout(35).Coef_A=1;
Coeffout(35).Coef_w=2*pi/365.25;
Coeffout(35).Coef_fi=pi/6;
Coeffout(35).Coef_Eps=0.1;

% (36:40) corresponding to In 2 with changing Seasonal phase of environment
% (for sensitivity analysis)
% Out 17 same to In 1
for i=1:3
Coeffout(i+35).name='out';
Coeffout(i+35).Coef_trend=1/365.25;
Coeffout(i+35).Coef_A=2;
Coeffout(i+35).Coef_w=2*pi/365.25;
Coeffout(i+35).Coef_fi=pi/6+pi/6*(i-1);
Coeffout(i+35).Coef_Eps=0.1;
end

for i=1:2
Coeffout(i+38).name='out';
Coeffout(i+38).Coef_trend=1/365.25;
Coeffout(i+38).Coef_A=2;
Coeffout(i+38).Coef_w=2*pi/365.25;
Coeffout(i+38).Coef_fi=pi/6-pi/6*i;
Coeffout(i+38).Coef_Eps=0.1;
end

Coeffout(41).name='out'; 
Coeffout(41).Coef_trend=1/365.25*0.2; % one fifth
Coeffout(41).Coef_A=2;
Coeffout(41).Coef_w=2*pi/365.25;
Coeffout(41).Coef_fi=pi/6;
Coeffout(41).Coef_Eps=0.1;

Coeffout(42).name='out';
Coeffout(42).Coef_trend=1/365.25;
Coeffout(42).Coef_A=0.4; % one fifth
Coeffout(42).Coef_w=2*pi/365.25;
Coeffout(42).Coef_fi=pi/6;
Coeffout(42).Coef_Eps=0.1;

%43 for 2 (make figure)
Coeffout(43).name='out';
Coeffout(43).Coef_trend=1/365.25*0.2; % one fifth
Coeffout(43).Coef_A=2*5; % fifth
Coeffout(43).Coef_w=2*pi/365.25;
Coeffout(43).Coef_fi=pi/6+pi/3; % additional pi/3
Coeffout(43).Coef_Eps=0.1;

%44 for 2 (make figure)
Coeffout(44).name='out';
Coeffout(44).Coef_trend=1/365.25*0.5; % half
Coeffout(44).Coef_A=2*2; % twice
Coeffout(44).Coef_w=2*pi/365.25;
Coeffout(44).Coef_fi=pi/6+pi/3; % additional pi/3
Coeffout(44).Coef_Eps=0.1;

%45 for 2 (make figure)
Coeffout(45).name='out';
Coeffout(45).Coef_trend=1/365.25*2; % twice
Coeffout(45).Coef_A=2*2; % twice
Coeffout(45).Coef_w=2*pi/365.25; 
Coeffout(45).Coef_fi=pi/6+pi/3; % additional pi/3
Coeffout(45).Coef_Eps=0.1;

%46 for 2 (make figure)
Coeffout(46).name='out';
Coeffout(46).Coef_trend=1/365.25*5; % one fifth
Coeffout(46).Coef_A=2*5; % twice
Coeffout(46).Coef_w=2*pi/365.25;
Coeffout(46).Coef_fi=pi/6+pi/3; % additional pi/3
Coeffout(46).Coef_Eps=0.1;

%47 for 8(in) (make figure)
Coeffout(47).name='out';
Coeffout(47).Coef_trend=1/365.25*0.2; % one fifth
Coeffout(47).Coef_A=2*5; % fifth
Coeffout(47).Coef_w=2*pi/365.25;
Coeffout(47).Coef_fi=pi/6+pi/3; % additional pi/3
Coeffout(47).Coef_Eps=1;

% 48 no signal of environment (even the noise)
Coeffout(48).name='out';
Coeffout(48).Coef_trend=0;
Coeffout(48).Coef_A=0;
Coeffout(48).Coef_w=0;
Coeffout(48).Coef_fi=0;
Coeffout(48).Coef_Eps=0;

save('Coeff','Coeffin','Coeffout')