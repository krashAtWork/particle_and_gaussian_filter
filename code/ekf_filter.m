% EECS 499 - Homework 3
% Question 1.2
% Author: Krashagi Gupta
% Case ID: kxg360
clear;
clc;

% initialise
%% 1.1
global alpha
global sigmaR
global sigmaPhi
sigmaR = 0.1;
sigmaPhi = 0.09;
alpha = [0.0001, 0.0001, 0.01, 0.0001, 0.0001, 0.0001 ];
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
for i = 1:6    
    plot(map(i,1),map(i,2),'r*')
    hold on ;
end

mu0 = [2; 2; 0]; 
zeta0 =  [10^-9 0 0; 0 10^-9 0;0 0 10^-9 ];
c0 = Zt(1,3);

mu_t = [0 ; 0 ; 0];
zeta_t = zeros(3,3);
temp_mu = mu0 ;
temp_zeta  =zeta0 ;
mu_collec = zeros(8, 3);
mu_bar_collec = zeros(8, 3);
zeta_collec = zeros(24, 3);
zeta_bar_collec = zeros(24, 3);

for k = 1: 8    
    [temp_mu ,temp_zeta, mu_t_bar,zeta_t_bar ] = ekf_algo(temp_mu, temp_zeta, Ut(k,:), Zt(k,:), Zt(k,3), map);
    mu_collec(k,:) = temp_mu' ;    
    zeta_collec(1 + 3*(k -1) : 3*k,:)  = temp_zeta ;
    mu_bar_collec(k,:) = mu_t_bar' ;    
    zeta_bar_collec(1 + 3*(k -1) : 3*k,:)  = zeta_t_bar ;
    
end
%  h1 = plot_gaussian_ellipsoid([1 1], [1 0.5; 0.5 1]);
%  h2 = plot_gaussian_ellipsoid([2 1.5], [1 -0.7; -0.7 1]);
%  h3 = plot_gaussian_ellipsoid([0 0], [1 0; 0 1]);
%  set(h2,'color','r'); 
%  set(h3,'color','g');

plot_everything(mu0', Ut, mu_collec, zeta_collec, mu_bar_collec, zeta_bar_collec )
% plot_gaussian_ellipsoid(m, C, sdwidth, npts, axh)
%% ekf function

function [mu_t ,zeta_t, mu_t_bar,zeta_t_bar] = ekf_algo(mu_t_1, zeta_t_1, ut, zt, ct, map)
%    test dimension
 

%    zt
%    ct
   map;
    global alpha
    global sigmaR
    global sigmaPhi
     del_t = 1;
     vt = ut(1,1);
     wt = ut(1,2); 
%     line 2
    theta = mu_t_1(3,1);
    
%    line 3
%    iniialisation     
    gt13 = (-vt/wt)* cos(theta) + (vt/wt)* cos(theta + (wt * del_t));
    gt23 =  (-vt/wt)* sin(theta) + (vt/wt)* sin(theta + (wt * del_t));
    Gt = [1 0 gt13 ; 0 1 gt23 ; 0 0 1]; %% 3* 3 matrix
    
%     line 4
%     initialisation
    vt11  = (-sin(theta)+ sin(theta + (wt*del_t)))/wt ;
    vt12 = (vt/wt^2)*(sin(theta)- sin(theta + (wt * del_t)))  +   (vt * del_t/wt)*(cos(theta + (wt * del_t)));
    vt21 = (cos(theta)- cos(theta + (wt*del_t)))/wt ;
    vt22 = (-vt/wt^2)*(cos(theta)- cos(theta + (wt * del_t)))  +   (vt * del_t/wt)*(sin(theta + (wt * del_t)));
    Vt = [vt11 vt12 ; vt21 vt22; 0 del_t]; %% 3*2 matrix
    
%     line 5
    Mt = [alpha(1)*vt^2 + alpha(2)*wt^2  0; 0 alpha(3)*vt^2 + alpha(4)*wt^2 ];% 2*2 matrix
    
%     line 6
    mu_t_bar = mu_t_1 + [gt23 ; -gt13 ; wt* del_t]; % 3*1 matrix
    
%     line 7
    zeta_t_bar = (Gt * zeta_t_1 * Gt') + (Vt * Mt * Vt');

    
%    line 8
    Qt = [sigmaR^2 0; 0 sigmaPhi^2]; %% 2*2 matrix
   
%     line 9 - > at one time instant only 1 observation
%       line 10
    j = ct;
%      line 11
      
    q = (map(j,1)- mu_t_bar(1,1))^2 + (map(j,2) - mu_t_bar(2,1))^2 ;% single value
    
%     line 12

    zt_cap = [sqrt(q);atan2(map(j,2) - mu_t_bar(2,1),map(j,1) - mu_t_bar(1,1)) - mu_t_bar(3,1)];% 2*1 %% looks wrong the theta value
%     zt

%     line 13
    Ht = [-(map(j,1)- mu_t_bar(1,1))/sqrt(q)   -(map(j,2) - mu_t_bar(2,1))/sqrt(q)  0;
            (map(j,2) - mu_t_bar(2,1))/q   -(map(j,1)- mu_t_bar(1,1))/q    -1 ];  % 2*3 matrix

%         line 14
    St = Ht * zeta_t_bar * Ht' + Qt ; % 2*2 

%     line 15
    %%% how to find inverse of matrix in matab
%     Y = inv(X) matrix must be square
    Kt = zeta_t_bar *   Ht' * inv(St) ;% 3*2
    
%     line 16 
   
%      zt   - zt_cap
    mu_t = mu_t_bar + Kt * (zt(:,1:2)' - zt_cap); % address concern % 3*3
    
%     line 17
%     how to make identity matrix of a given size in matlab
%     I = eye(n)
    zeta_t = (eye(3) - Kt*Ht)* zeta_t_bar;
end
%% plot everything
function [] = plot_everything(x0,Ut, mu_collec,zeta_collec,mu_bar_collec, zeta_bar_collec)   
x0
Ut
mu_collec
zeta_collec
mu_bar_collec
zeta_bar_collec
%     x = 2;
%     y = 2;
%     theta = 0 ;
%     v = pi/2;
%     w = pi/2;
    X = zeros(8, 3)
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

%     for k = 1 : 1
%          plot_gaussian_ellipsoid( mu_bar_collec(k,:) , zeta_bar_collec(1 + 3*(k -1) : 3 * k,:))        
%     end
    %  h1 = plot_gaussian_ellipsoid([1 1], [1 0.5; 0.5 1]);
%  h2 = plot_gaussian_ellipsoid([2 1.5], [1 -0.7; -0.7 1]);
%  h3 = plot_gaussian_ellipsoid([0 0], [1 0; 0 1]);
%  set(h3,'color','g');
    
    for k = 1 : 8 
        plot_gaussian_ellipsoid( mu_collec(k,:) , zeta_collec(1 + 3*(k -1) : 3 * k,:));
    end
  
end
%% functions
function [new_posi] = move(old_posi, miu, delta_t)
    new_posi = old_posi + [-miu(1)/miu(2) * sin(old_posi(3)) + miu(1)/miu(2) * sin(old_posi(3) + miu(2) * delta_t), miu(1)/miu(2) * cos(old_posi(3)) - miu(1)/miu(2) * cos(old_posi(3) + miu(2) * delta_t), miu(2) * delta_t];
end