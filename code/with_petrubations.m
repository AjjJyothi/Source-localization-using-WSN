clc;
close all;
clear all;% N sensors, 1 source. Using sensor 1 as reference i.e (x1=0, y1=0, z1=0)
% -----------------------------------------------------
% Definition
% -----------------------------------------------------
nRun=100; % number of Monte Carlo runs
% uncomment one of them
 bML=1; % turn off ML calculation
%bML=1; % turn on ML calculation
% uncomment one of them
perturb=1; % turn off location perturbation
% perturb=1; % turn on location perturbation
% ----------------------------------------------------------------
% Actual source location (m) in Cartesian coordinates x, y and z
% Note: For simplicity, we only varies y for our simulation
% ----------------------------------------------------------------
xs_src_actual=[0];
% Varies the Y position (Choose 1)
%------------------------------------
ys_src_actual=[100]; %100 m
zs_src_actual=[0];
Rs_actual=sqrt(xs_src_actual.^2 + ys_src_actual.^2 + zs_src_actual.^2);
% calculate corresponding range Rs
bearing_actual=[xs_src_actual; ys_src_actual; zs_src_actual]/Rs_actual;
% calculate corresponding bearing
% ----------------------------------------------------------------
% Actual sensor location (m) in Cartesian coordinates x, y and z
% Note: For simplicity, we only use integers and then multiply with
% a scaling factor to produce the actual coordinates.
% e.g. [5 10 15] = [1 2 3] * 5 ;
% ----------------------------------------------------------------
% Scale wrt to 1m (Choose 1)
%------------------------------
scale_dist = 10; % 10 m
% (Choose 1 of the following sensor configuration for study)
%-------------------------------------------------------------
% 12x2 sensors arranged 2 rows
  xi=[0 0 1 1 2 2 3 3 4 4 5 5 0 0 1 1 2 2 3 3 4 4 5 5 ].*scale_dist;
  yi=[0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1].*scale_dist;
  zi=[0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1].*1.0;
 %20x2 sensors arranged 2 rows
%  xi=[0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9].*scale_dist;
%  yi=[0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1].*scale_dist;
%  zi=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1].*1.0;
% Soldier configuration (each with 2 sensors). z= 1m apart vertically
temp=size(xi);
nSen=temp(1,2); % number of sensor (>4)
noisestd=1;
if (perturb==1)
randn('state',0);
tmp1=randn(3, nSen);
for i=1:2:nSen
xi(i)= xi(i) + noisestd*tmp1(1,i);
xi(i+1)= xi(i+1)+(noisestd+0.01)*tmp1(1,i); % less variance on the body
yi(i)= yi(i) + noisestd*tmp1(2,i);
yi(i+1)= yi(i+1)+(noisestd+0.01)*tmp1(2,i); % less variance on the body
end
zi=zi+0.01*tmp1(3,:); % less variance on the body
end
% RD noise (Choose 1)
% -----------------------------------------------------
Noise_Factor=0.02;
 % noise std = Std_Norm * (source distance). %we expect bigger noise variance for larger distance.
Noise_Var=(Noise_Factor*Rs_actual)^2;
% -----------------------------------------------------
% Functions
% -----------------------------------------------------
% Random Process
% AWGN
randn('state',0);
noise = sqrt(Noise_Var)*randn(nRun, 1);
%noise_mean = mean(noise, 2); % average along row
for k=1:nRun % Monte Carlo Simulation
Xi=[xi' yi' zi'];
Di= sqrt ((xi-xs_src_actual).^2 + (yi-ys_src_actual).^2 + (zi-zs_src_actual).^2);%distance between source ans sensor i
Ri= sqrt ((xi).^2 + (yi).^2 + (zi).^2);%distance between origin and sensor i
locSen=[xi' yi' zi'];
% using N sensors
for i=1:nSen-1
%d21=Di(2)-Di(1);
%d31=Di(3)-Di(1);..
%dn1=Di(n)-Di(1);
%d=[d21;d31;...;dn1];
d(i,1)=Di(i+1)-Di(1)+noise(k); %add noise to RD estimates%distance between sensor i,1
% delta2=Ri(2)^2-d(1)^2;
% delta3=Ri(3)^2-d(1)^2;...
% deltan=Ri(n)^2-d(1)^2;
% delta=[delta2;delta3;...deltan];
delta(i,1)=Ri(i+1)^2-d(1)^2;%delta
% s2= [xi(2) yi(2) zi(2)];
% s3= [xi(3) yi(3) zi(3)];...
% sn= [xi(n) yi(n) zi(n)];
% s=[s2;s3;...sn];
s(i,:)=[xi(i+1) yi(i+1) zi(i+1)];
end
% define weight (positive definite, i.e., diagonal positive and symmetrical)
w=eye(nSen-1); % set to identity matrix for unweighted case
Sw=(s'*w*s)^(-1)*s'*w;
Ps=s*Sw;
Ps_ortho=eye(nSen-1)-Ps;%projection matrix
%------------------------------------------------
%SI method
%------------------------------------------------

Rs_SI_cal=0.5*(d'*Ps_ortho*w*Ps_ortho*delta)/(d'*Ps_ortho*w*Ps_ortho*d);%s cap
% Calculate Xs for SI method
Xs_row_SI = 0.5*Sw*(delta-2*Rs_SI_cal*d);
Xs_SI(k,:)=Xs_row_SI' ;
Rs_SI(k,:)=sqrt(Xs_SI(k,1)^2 + Xs_SI(k,2)^2 + Xs_SI(k,3)^2);%range estimate,norm(s cap)
bearing_SI(k,:)=Xs_SI(k,:)/Rs_SI(k,:);
% error_row= delta - 2*Rs_SI*d - 2*s*Xs_row_SI; %error
% error(k,:)=error_row';
% %------------------------------------------------------
% % Maximum Likelihood Method
% % Objective function is contained in mlobjfun.m
% %----------------------------------------------------------
if (bML==1)
x0 = Xs_SI(k,:); % As value obtained from SI as starting guess
 %x0 = [0 ys_src_actual 0]; % Starting guess
options = optimset('LargeScale','off');
% % LevenbergMarquardt
%options=optimset(options,levenbergMarquardt','on'); % LM
 %options=optimset(options, 'levenbergMarquardt','off'); % Gauss Newton
 options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');
 [x,resnorm,residual,exitflag,output]=lsqnonlin(@(x) mlobjfun(x,locSen,Noise_Var,d),x0,[],[],options);
Xs_ML(k,:)=x;
Rs_ML(k,:)=sqrt(Xs_ML(k,1)^2+Xs_ML(k,2)^2+Xs_ML(k,3)^2);
bearing_ML(k,:)=Xs_ML(k,:)/Rs_ML(k,:);
end

% ---------------------------------------------------------------------
% Calculate bias (i.e., errors) for source location, range and bearing
% estimates
% ---------------------------------------------------------------------
% Calculate mean
% -----------------------------------------------------
% SI
bias_Xs_SI(k,1)=Xs_SI(k,1)-xs_src_actual;
bias_Xs_SI(k,2)=Xs_SI(k,2)-ys_src_actual;
bias_Xs_SI(k,3)=Xs_SI(k,3)-zs_src_actual;
% ML
if (bML==1)
bias_Xs_ML(k,1)=Xs_ML(k,1)-xs_src_actual;
bias_Xs_ML(k,2)=Xs_ML(k,2)-ys_src_actual;
bias_Xs_ML(k,3)=Xs_ML(k,3)-zs_src_actual;
end
end
bias_Rs_SI = Rs_SI-Rs_actual;
bias_bearing_SI = 180/pi*acos(bearing_SI*bearing_actual);
if (bML==1)
bias_Rs_ML=Rs_ML-Rs_actual;
bias_bearing_ML = 180/pi*acos(bearing_ML*bearing_actual);
end
meanxs_SI=mean(bias_Xs_SI(:,1));
meanys_SI=mean(bias_Xs_SI(:,2));
meanzs_SI=mean(bias_Xs_SI(:,3));
meanrs_SI=mean(bias_Rs_SI);
meanbear_SI=mean(bias_bearing_SI);
vect_mean_SI=[meanxs_SI;meanys_SI;meanzs_SI;meanrs_SI;meanbear_SI];
%ML
if (bML==1)
meanxs_ML=mean(bias_Xs_ML(:,1));
meanys_ML=mean(bias_Xs_ML(:,2));
meanzs_ML=mean(bias_Xs_ML(:,3));
meanrs_ML=mean(bias_Rs_ML);
meanbear_ML=mean(bias_bearing_ML);
vect_mean_ML=[meanxs_ML;meanys_ML;meanzs_ML;meanrs_ML;meanbear_ML];
end
% Calculate Variance = E[(a - mean)^2]
% -----------------------------------------------------
varxs_SI=var(bias_Xs_SI(:,1));
varys_SI=var(bias_Xs_SI(:,2));
varzs_SI=var(bias_Xs_SI(:,3));
varrs_SI=var(bias_Rs_SI);
varbear_SI=var(bias_bearing_SI);
vect_var_SI=[varxs_SI;varys_SI;varzs_SI;varrs_SI;varbear_SI];
%ML
if (bML==1)
varxs_ML=var(bias_Xs_ML(:,1));
varys_ML=var(bias_Xs_ML(:,2));
varzs_ML=var(bias_Xs_ML(:,3));
varrs_ML=var(bias_Rs_ML);
varbear_ML=var(bias_bearing_ML);
vect_var_ML=[varxs_ML;varys_ML;varzs_ML;varrs_ML;varbear_ML];
end
% Calculate second moment (RMS)= sqrt {E[a^2]} = sqrt {mean^2 + variance}
% -----------------------------------------------------
rmsxs_SI=sqrt(mean(bias_Xs_SI(:,1)).^2+varxs_SI);
rmsys_SI=sqrt(mean(bias_Xs_SI(:,2)).^2+varys_SI);
rmszs_SI=sqrt(mean(bias_Xs_SI(:,3)).^2+varzs_SI);
rmsrs_SI=sqrt(mean(bias_Rs_SI).^2+varrs_SI);
rmsbear_SI=sqrt(mean(bias_bearing_SI).^2+varbear_SI);
vect_rms_SI=[rmsxs_SI;rmsys_SI;rmszs_SI;rmsrs_SI;rmsbear_SI];
%ML
if (bML==1)
rmsxs_ML=sqrt(mean(bias_Xs_ML(:,1)).^2+varxs_ML);
rmsys_ML=sqrt(mean(bias_Xs_ML(:,2)).^2+varys_ML);
rmszs_ML=sqrt(mean(bias_Xs_ML(:,3)).^2+varzs_ML);
rmsrs_ML=sqrt(mean(bias_Rs_ML).^2+varrs_ML);
rmsbear_ML=sqrt(mean(bias_bearing_ML).^2+varbear_ML);
vect_rms_ML=[rmsxs_ML;rmsys_ML;rmszs_ML;rmsrs_ML;rmsbear_ML];
end
% Calculate Cramer Rao Bound
%
cov_mat=Noise_Var.*(0.5*ones(length(d))+0.5*eye(length(d)));
for i=1:length(d)
a1=[xs_src_actual-locSen(i+1,1) ys_src_actual-locSen(i+1,2) zs_src_actual-locSen(i+1,3)];
a2=sqrt((xs_src_actual-locSen(i+1,1))^2+(ys_src_actual-locSen(i+1,2))^2+(zs_src_actual-locSen(i+1,3))^2);
b1=[xs_src_actual-locSen(1,1) ys_src_actual-locSen(1,2) zs_src_actual-locSen(1,3)];
b2=sqrt((xs_src_actual-locSen(1,1))^2+(ys_src_actual-locSen(1,2))^2+(zs_src_actual-locSen(1,3))^2);
jacobian(i,:)= (a1/a2)-(b1/b2);
end
fisher=jacobian'*inv(cov_mat)*jacobian;
crlb= trace(fisher^-1); % compare with MSE of Rs
% -----------------------------------------------------
% Generate Plots
% -----------------------------------------------------
% hfig1=figure;
if (bML==1)
plot(xi, yi,'kv', Xs_SI(:,1), Xs_SI(:,2),'kx',Xs_ML(:,1), Xs_ML(:,2), 'ko', xs_src_actual, ys_src_actual, 'r^'); % plot both SI and ML
legend('sensor location', 'calculated source location(SI)','calculated source location (ML)', 'actual source location ');
else
plot(xi, yi,'kv', Xs_SI(:,1), Xs_SI(:,2), 'kx', xs_src_actual, ys_src_actual, 'r^');
legend('sensor location', 'calculated source location (SI)', 'actual source location');
% plot just SI only
end
title('Sensor and Source Location');
xlabel('Distance (metres) in X direction');
ylabel('Distance (metres) in Y direction');
% axis([-10 30 0 120]);
% variance=Noise_Var
% scaling=scale_dist
% position=ys_src_actual
 rms=[rmsxs_SI, rmsys_SI, rmszs_SI,rmsrs_SI, rmsbear_SI]
 mean=[meanxs_SI, meanys_SI, meanzs_SI,meanrs_SI, meanbear_SI]
 variance=[varxs_SI, varys_SI, varzs_SI,varrs_SI, varbear_SI]
  if (bML==1)
  rms=[rmsxs_ML, rmsys_ML, rmszs_ML,rmsrs_ML, rmsbear_ML]
  mean=[ meanxs_ML, meanys_ML,meanzs_ML, meanrs_ML, meanbear_ML]
  variance=[varxs_ML, varys_ML,varzs_ML, varrs_ML, varbear_ML]
 end

