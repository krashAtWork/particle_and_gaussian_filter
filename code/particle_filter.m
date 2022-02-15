% EECS 499 - Homework 3
% Question 1.2
% Author: Krashagi Gupta
% Case ID: kxg360
clear;
clc;

% initialise
%% 1.2
global M
M = 100;
x0 = [2, 2, 0]; % vector inital pose
X0 = repmat(x0, M, 1);

Ut = [[1,0.0001]; [1,0.0001];[1,0.0001]; [pi/2, pi/2]; [pi/2, pi/2];[1,0.0001]; [1,0.0001];[1,0.0001]];

Zt = [[2.276, 5.249, 2];[4.321, 5.834, 3];[3.418, 5.869, 3];[3.774, 5.911, 4];
    [2.631, 5.140, 5];[4.770,5.791,6];[3.828, 5.742, 6];[3.153, 5.739, 6]];

m1 = [0, 0];
m2 = [4, 0];
m3 = [8, 0];
m4 = [8, 6];
m5 = [4, 6];
m6 = [0, 6];
map = [m1; m2; m3; m4; m5; m6];
h = figure;
for i = 1:6    
    plot(map(i,1),map(i,2),'r*')
    hold on ;
end
Xr = zeros(8 * M,4);
Xbar = zeros(8 * M,4);
step = 0;
tempX = zeros(M,3);
tempX = X0;

for q = 1 : 8
     [Xr(1+ M*(q - 1):q*M,:),Xbar(1+ M*(q - 1):q*M,:)]= mcl_algo(tempX, Ut(q,:), Zt(q,:), map);  %%%% problem is here for 2nd one
     tempX  = Xr(1+ M*(q - 1):q * M,:);    
end

plot_everything(x0, Ut ,Xr, Xbar)

%% mcl algorithm
function [Xt, Xt_bar] = mcl_algo(Xt_1, ut, zt, map)
    
   global M
   alpha = [0.0001, 0.0001, 0.01, 0.0001, 0.0001, 0.0001 ];

    % line 2
    Xt_bar =  zeros(M,4);
    Xt = zeros(M,4);
    tempX = zeros(1,3);
    
    % line 3
    for m = 1: M
   

     tempX = sample_motion_model_velocity1(ut, Xt_1(m,:))';
           
%            
%    line 5
       tempWt = measurement_model(zt, tempX, map); %%%% have to solve
%    line 6

        Xt_bar(m,:)= [tempX tempWt] ;
    end
     % rearrange by wt
     
         Xt_bar = sortrows(Xt_bar,1);
%         Xt_bar(m,:)
%       line 7
        Xt = resample(Xt_bar);
        
   
   
end
%% sample_motion_model

function [xn]= sample_motion_model(ut0, xp, alpha)

    
        xn =[0,0,0];
        xo = xp(1);
        yo = xp(2);
        thetao = xp(3);
        delT = 1 ;

    %     building an error array
        errorArray = sample_error(ut0,alpha);
    
        vcp  = ut0(1) + errorArray(1,1);
        wcp = ut0(2) + errorArray(2,1);
        ycp = errorArray(3,1);


        xn(1,1) = xo  + (-vcp/wcp)* sin(thetao) + (vcp/wcp)*sin(thetao + (wcp* delT)); 
        xn(1,2) = yo  + (vcp/wcp)* cos(thetao) - (vcp/wcp)*cos(thetao + (wcp* delT));
        xn(1,3) = thetao  +  wcp*delT  + ycp*delT ;
      
%     end
end

%% sample error 

function[errorArray] = sample_error(ut, alp)

errorArray = [0; 0; 0;];
v = ut(1);
w = ut(2);
vcpE = sample(v, w, alp(1), alp(2));
wcpE = sample(v, w, alp(3), alp(4));
ycpE =  sample(v, w, alp(5), alp(6));
errorArray(1,1)= vcpE;
errorArray(2,1)= wcpE;
errorArray(3,1)= ycpE;

end
%% sample
function[err] = sample(v, w, alfF, alfS)
varnc = alfF*v + alfS*w;
mu = 0;
sigma = sqrt(varnc);
err = normrnd(mu,sigma);
end
%% sample2
function [xt]=sample_motion_model_velocity1(ut,xtm1)
N = 1;
alpha=[0.0001;0.0001;0.01;0.0001;0.0001;0.0001];
vv=ut(1)+normrnd(0,(alpha(1)*ut(1)^2+alpha(2)*ut(2)^2)^0.5,[N,1]);
ww=ut(2)+normrnd(0,(alpha(3)*ut(1)^2+alpha(4)*ut(2)^2)^0.5,[N,1]);
rr=normrnd(0,(alpha(5)*ut(1)^2+alpha(6)*ut(2)^2)^0.5,[N,1]);
    for i=1:N
    xt(1,:)=xtm1(1)-vv./ww.*sin(xtm1(3))+vv./ww.*sin(xtm1(3)+ww.*1);
    xt(2,:)=xtm1(2)+vv./ww.*cos(xtm1(3))-vv./ww.*cos(xtm1(3)+ww.*1);
    xt(3,:)=xtm1(3)+ww.*1+rr.*1;
    end
end


%%
function [reqProb] = measurement_model(z,X,m)
      % input received right
%     input form 
%     z -> [d,a]
      d = z(1);
      a = z(2);
      s = z(3);
%     X -> [x, y, theta]
      x = X(1);
      y = X(2);
      theta = X(3);
%     m -> [mx, my]
      mx = m(s,1);
      my = m(s,2);
      % noise parameters
      sigmaR = 0.1;
      sigmaPhi = 0.09;
  
    dCp = sqrt((mx - x)^2 + (my - y)^2);
    alphaCp =  atan2(my - y,mx - x) - theta;
%     pDet = prob(dCp - d, sigmaR^2)* prob(alphaCp - a, sigmaPhi^2);
    pDet = prob(dCp - d, sigmaR^2)* prob(alphaCp - a, sigmaPhi^2);
  
    Zdet = 1;
    reqProb = Zdet * pDet;
%     1/(2*pi*sigmar^2)^0.5* exp(-0.5*(zt(1)-rhat)^2/sigmar^2)
    
    
end
%% prob function
function [ansr] = prob(a, var)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
ansr =  (1/(2*pi*var)^0.5)* exp(-0.5*(a)^2/var);
end

%% plot nominal motion
function [] = plot_everything(x0,Ut, Xr, Xbar)   
%     x = 2;
%     y = 2;
%     theta = 0 ;
%     v = pi/2;
%     w = pi/2;
    X = zeros(8, 3);
    global M;
    global delta_t;
    delta_t = 1;
    old_location = x0;
    for t = 1: 8        
        temp_X = move(old_location, Ut(t,:), delta_t);
        line([old_location(1), temp_X(1)], [old_location(2), temp_X(2)], 'Color','black');
        hold on ;
       
        X(t,:) = temp_X;
        old_location = temp_X;
    end
    
    for t = 1:8
        plot(X(t,1),X(t,2),'r*'); 
        hold on ;
    end   
    for g = 1 : 8 * M
        scatter(Xbar(g,1),Xr(g,2),10,'b','filled');
        hold on;
    end
    for g = 1 : 8 * M
        scatter(Xr(g,1),Xr(g,2),10,'k','filled');
        hold on;
    end 
    

end

%% functions
function [new_posi] = move(old_posi, miu, delta_t)
    new_posi = old_posi + [-miu(1)/miu(2) * sin(old_posi(3)) + miu(1)/miu(2) * sin(old_posi(3) + miu(2) * delta_t), miu(1)/miu(2) * cos(old_posi(3)) - miu(1)/miu(2) * cos(old_posi(3) + miu(2) * delta_t), miu(2) * delta_t];
end

%% resample
function [Xt] = resample(Xt_bar)
global M
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    dim = size(Xt_bar); 
    Xt = zeros(dim); % ensuring output has same dimensions as input
    Xt_bar_norm =  Xt_bar./sum(Xt_bar);
    cumSum =  cumsum(Xt_bar_norm);
    for j = 1: M         
         chance = rand(1, 1);
         % works right
         for k = 1 : M             
            if(cumSum(k,1) >= chance)                
                 Xt(j,:) = Xt_bar(k,:);
                 break;
            end        
         end
    end
end

 %% shhd  
function [q]=landmark_model_known_correspondence(zt,ct,xt,m)
j=ct;
rhat=((m(1,j)-xt(1))^2+(m(2,j)-xt(2))^2)^0.5;
phihat=atan2(m(2,j)-xt(2),m(1,j)-xt(1))-xt(3);
q=1/(2*pi*0.01)^0.5*exp(-0.5*(zt(1)-rhat)^2/0.01)*1/(2*pi*0.09^2)^0.5*exp(-0.5*(zt(2)-phihat)^2/0.09^2);

end


