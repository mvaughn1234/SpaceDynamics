%% Mathew Vaughn - PS8 EKF Code
% Template Courtesy of Dr. Becnel
function [yHatArray,PArray,debugStruct] = EKF(tVec_sec,...
    R,Q,yHatVec_t0,P_t0,xArray,param)

% Run parameters:
nMeas = length(tVec_sec);
ny = size(yHatVec_t0,1);
nx = size(xArray,1);
% Initial Storage:
yHatArray = zeros(ny,nMeas);
PArray = zeros(ny,ny,nMeas);
yHatArray(:,1) = yHatVec_t0;
PArray(:,:,1) = P_t0;
%
debugStruct.innovTestStat = zeros(1,nMeas);
debugStruct.yBarVec_tkp1 = zeros(ny,nMeas);
debugStruct.HArray = zeros(nx,ny,nMeas);
debugStruct.PhiArray = zeros(ny,ny,nMeas);
debugStruct.hVecMat = zeros(nx,nMeas);
debugStruct.KArray = zeros(ny,nx,nMeas);
% Initialize for loop:
yHatVec_tk = yHatVec_t0;
P_tk = P_t0;
% Loop through all measurements:
for k = 1:(nMeas-1) %No measurement at time 0.
   % Time: tk to tkp1
   tk_sec = tVec_sec(k);
   tkp1_sec = tVec_sec(k+1);
   %
   % Dynamics Model:
   %    This should take the current best state estimate (yHat at tk) and
   %    propagate it to the next time step t(k+1).  This propagated state
   %    is termed yBar at time t(k+1) and is not the final state estimate
   %    for time-step t(k+1) because it has not yet had the measurement
   %    associated with this time processed to make a more optimal
   %    estimate.  It will also produce Phi(t(k),t(k+1)).  
   %    The parameter structure should contain fields that are used by the
   %    dynamics model (for example, muE:   param.muE = 398600).
   %    
   [yBarVec_tkp1,Phi_k2kp1] = dynamicsModel(yHatVec_tk,...
       tk_sec,tkp1_sec,param);
   %
   % Measurement Model:
   %    This should take the best state estimate that has been propagated
   %    from tk to tk+1 and produce the measurements expected from that
   %    state.  It will also take in the parameters necessary to evaluate
   %    the measurement model.
   %
   [hVec, H] = measurementModel(yBarVec_tkp1,tkp1_sec,param);
   %
   %
   % Estimation Equations:
   %
   % Propagate the covariance to t(k+1) using Phi and Q:
   PBar_tkp1 = Phi_k2kp1*P_tk*Phi_k2kp1'+Q;
   % Compute the Kalman gain using PBar (above), H, and R:
   K = (PBar_tkp1*H')*(H*PBar_tkp1*H'+R)^-1;
   % Compute the difference between the actual measurement and the expected 
   % measurement at time t(k+1).  Use hVec, not H*y.
   errorTermVec = xArray(:,k+1) - hVec;
   % Compute the state estimate at time t(k+1) using the propagated state
   % estimate from time t(k) to t(k+1) and then the correction term of the
   % Kalman gain times the errorTermVec:
   yHatVec_tkp1 = yBarVec_tkp1 + K*errorTermVec;
   % Compute the new covariance estimate using the propagated covariance
   % Pbar, and H,R,K:
   P_tkp1 = (eye(ny)-K*H)*PBar_tkp1;
   %
   %   
   % Setup the next iteration: (tkp1->tk)
   yHatArray(:,k+1) = yHatVec_tkp1;
   PArray(:,:,k+1) = P_tkp1;
   yHatVec_tk = yHatVec_tkp1;
   P_tk = P_tkp1;
   %
   % Store data for debugging:
%    debugStruct.innovTestStat(k+1) =  ???;
   debugStruct.yBarVec_tkp1(:,k+1) = yBarVec_tkp1;
   debugStruct.HArray(:,:,k+1) = H;
   debugStruct.PhiArray(:,:,k+1) = Phi_k2kp1;
   debugStruct.hVecMat(:,k+1) = hVec;
   debugStruct.KArray(:,:,k+1) = K;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dynamicsModel:
%
% This function will take the state estimate and propagate it from t(k) to 
% t(k+1) and produce Phi(t(k),t(k+1)).  The distinction between
% linear and nonlinear problems is important, assuming no control or
% process noise:
%   Linear:     yBarVec at t(k+1) = Phi(tk,tk+1) * yHatVec at t(k)
%   Nonlinear:  yBarVec at t(k+1) = f[yHatVec at t(k), tk, tk+1]  In this
%       case the state is NOT Phi(tk,tk+1) * yHatVec at t(k).  Instead, 
%       you should use the state out of the propagator that solved the ode.
%       In the special case of linear dynamics:
%           f[yHatVec at t(k), tk, tk+1] = Phi(tk,tk+1) * yHatVec at t(k)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yBarVec_tkp1,Phi_k2kp1] = dynamicsModel(yVec_tk,...
    tk_sec,tkp1_sec,param)
t_range = [tk_sec:tkp1_sec];
[t2,y,stm] = Position_2BP_STM(yVec_tk,t_range,10E-13,param.Re,0,param.muE);
yBarVec_tkp1 = y(end,:)';
Phi_k2kp1 = squeeze(stm(end,:,:));
Phi_k2kp1 = reshape(Phi_k2kp1,[6 6]);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% measurementModel:
%
% This function will take the propagated state (yBarVec at t(k+1)) and
% compute the expected measurement at that time.  This measurement is
% typically referred to as xBarVec at t(k+1) or hVec.
% The distinction between linear and nonlinear problems:
%   Linear:     hVec = H*yBarVec at t(k+1)
%   Nonlinear:  hVec = h[yBarVec at t(k+1), tk+1]  In this
%       case the measurement is NOT: H * yBarVec at t(k+1).  
%       Instead, you should use the actual measurement as hVec.
%       In the special case of linear dynamics:
%           h[yBarVec at t(k+1), tk+1] = H*yBarVec at t(k+1)
%
% Note: the partial derivatives can be computed using MATLAB's symbolic
% toolbox or by hand.  It is recommended to check the Jacobian numerically.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hVec, H] = measurementModel(yBarVec_tkp1,tkp1,param)
arg = ((2*pi)/(60*60*24))*tkp1;
Rz = [cos(arg) -sin(arg) 0;sin(arg) cos(arg) 0;0 0 1];
stat1_rot = Rz*param.stat1;
stat2_rot = Rz*param.stat2;
stat3_rot = Rz*param.stat3;
stat4_rot = Rz*param.stat4;

stat1_meas = norm(yBarVec_tkp1(1:3)-stat1_rot);
stat2_meas = norm(yBarVec_tkp1(1:3)-stat2_rot);
stat3_meas = norm(yBarVec_tkp1(1:3)-stat3_rot);
stat4_meas = norm(yBarVec_tkp1(1:3)-stat4_rot);

H = [(yBarVec_tkp1(1)-stat1_rot(1))/stat1_meas (yBarVec_tkp1(2)-stat1_rot(2))/stat1_meas (yBarVec_tkp1(3)-stat1_rot(3))/stat1_meas 0 0 0;
    (yBarVec_tkp1(1)-stat2_rot(1))/stat2_meas (yBarVec_tkp1(2)-stat2_rot(2))/stat2_meas (yBarVec_tkp1(3)-stat2_rot(3))/stat2_meas 0 0 0;
    (yBarVec_tkp1(1)-stat3_rot(1))/stat3_meas (yBarVec_tkp1(2)-stat3_rot(2))/stat3_meas (yBarVec_tkp1(3)-stat3_rot(3))/stat3_meas 0 0 0;
    (yBarVec_tkp1(1)-stat4_rot(1))/stat4_meas (yBarVec_tkp1(2)-stat4_rot(2))/stat4_meas (yBarVec_tkp1(3)-stat4_rot(3))/stat4_meas 0 0 0];

hVec = [stat1_meas stat2_meas stat3_meas stat4_meas]';

return