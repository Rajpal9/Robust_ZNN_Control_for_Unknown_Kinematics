%% This is a program for unknown jacobian with constraints on control and performance constraints
clear all

addpath("Kinematics\")
addpath("Paths\")


%% Basic Parameters
% time parameters
dt = 0.001; % time step
t_end = 10;% end time
t_sim = t_end; %simulation time for graphs
t_span = 0:dt:t_end; %time vector for simulation

% robot prelims
m = 3; % no. of equations
n = 7; % no. of joints
% robot params
%% parameters of kuka iiwa robot
% uncertian parameters
a8_unc = 0.126; %length of the end effector
d1_unc = 0.36;
d3_unc = 0.42;
d5_unc = 0.4;




%% Initial Conditions
% state
th(:,1) = [-0.2578,0.2946,0.5946,-0.9896,0.4993,0.3100,0.8620];
th_unres(:,1) = th(:,1);
%th(:,1) = [pi/4;3*pi/4;pi/3;pi/10;pi/2];
th_dot(:,1) = [0;0;0;0;0;0;0];




%% ZNN parameters
% observer ZNN gain
gamma_o = 100;
%controller ZNN gain
gamma_c = 100;


%% constraints


% for cardiod
th_dot_pos = 2*[0.1;0.1;0.1;0.1;0.1;0.1;0.1];
th_dot_neg = -th_dot_pos;
th_pos = 1.5*ones(n,1);
th_neg = -th_pos;
%%for paper
delta_con = 100;
p_0_con = 0.2;
p_inf_con = 5e-3;

delta_obs = 100;
p_0_obs = 0.4;
p_inf_obs = 5e-3;
k_pos = 5;
k_int = 00; %
kappa = 1;




% for cool_2
% th_dot_pos = 2*[0.1;0.1;0.1;0.1;0.1;0.1;0.1];
% th_dot_neg = -th_dot_pos;
% th_pos = 1.5*ones(n,1);
% th_neg = -th_pos;
% kappa = 3;

%th_pos = deg2rad([170;120;170;120;170;120;175]); % actual joint limits
%th_neg = - th_pos;

% th_pos =  [-0.17;1;0.65;-0.86;0.80;0.32;0.2]
% th_neg = [-0.38;1;1;-1.2;-0.5;0.11;-0.2]
th_dot_pos_plot = th_dot_pos*ones(1,size(t_span,2)); %upper bound for plot
th_dot_neg_plot = th_dot_neg*ones(1,size(t_span,2)); % lower bound for plot;

th_pos_plot = th_pos*ones(1,size(t_span,2)); %upper bound for plot
th_neg_plot = th_neg*ones(1,size(t_span,2)); % lower bound for plot;



for i  = 1:length(t_span)
    t = t_span(i);
    % magnitude of uncertianities
    unc_a8 = 0.01;
    unc_d1 = 0.03 + 0.07*cos(0.1*t);
    unc_d3 = 0.03 + 0.06*cos(0.1*t);
    unc_d5 = 0.04;

    % actual parameters
    a8 = a8_unc + unc_a8; %length of the end effector
    d1 = d1_unc + unc_d1;
    d3 = d3_unc + unc_d3;
    d5 = d5_unc + unc_d5;
    if i==1
        r_actual(:,1) = forward_map_kuka(a8,d1,d3,d5,th(:,1));
        dev = 0*[0.01;0.01;0.01];
        r_desired(:,1) = r_actual(:,1)+dev;

        r_aint(:,1) = r_actual(:,1);
        r_dint(:,1) = r_desired(:,1);
        
    end
    J((i-1)*m+1:i*m,:) = Jacobian_matrix_kuka(a8,d1,d3,d5,th(:,i));
    J_dot((i-1)*m+1:i*m,:) = Jacobian_dot_kuka(a8,d1,d3,d5,th(:,i),th_dot(:,i));
    
    %% the path
    % paper sim
     c = 0.029;
     tilt =  pi/4;
     shape = 'cardioid';

    [r_desired(:,i),rd_dot(:,i),rd_ddot(:,i)] = path_3D(c,tilt,r_desired(:,1),t,t_end,shape);
     r_actual(:,i) = forward_map_kuka(a8,d1,d3,d5,th(:,i));
  

    for k  = 1:n
        if th_dot_pos(k) <= kappa*(th_pos(k)-th(k,i))
            nu_pos(k,i) = th_dot_pos(k);
            nu_pos_dot(k,i) = 0;
        else
            nu_pos(k,i) = kappa*(th_pos(k)-th(k,i));
            nu_pos_dot(k,i) = -kappa*th_dot(k,i);           
        end
        
        if th_dot_neg(k) >= kappa*(th_neg(k)-th(k,i))
            nu_neg(k,i) = th_dot_neg(k);
            nu_neg_dot(k,i) = 0;
        else
            nu_neg(k,i) = kappa*(th_neg(k)-th(k,i));
            nu_neg_dot(k,i) = -kappa*th_dot(k,i);           
        end
    end
    
    del = diag(nu_pos(:,i)-nu_neg(:,i));
    del_dot = diag(nu_pos_dot(:,i)-nu_neg_dot(:,i));
    
    
    if i == 1
        xi_u(:,i) = log((nu_pos(:,i)-th_dot(:,i))./(th_dot(:,i)-nu_neg(:,i)));
        g_xi(:,i) = exp(xi_u(:,i))./(1+exp(xi_u(:,i)));
       J_hat((i-1)*m+1:i*m,:) = Jacobian_matrix_kuka(a8_unc,d1_unc,d3_unc,d5_unc,th(:,i))
       J_hat_dot = round(J_dot((i-1)*m+1:i*m,:),1);
        J((i-1)*m+1:i*m,:),J_dot((i-1)*m+1:i*m,:)
        r_aint(:,i) = r_actual(:,i);
        r_dint(:,i) = r_desired(:,1);
    else
        r_aint(:,i) = r_actual(:,i-1)*dt + r_aint(:,i-1);
        r_dint(:,i) = r_desired(:,i-1)*dt + + r_dint(:,i-1);
    end
    
    % performance criteria
    p_con(i) = (p_0_con-p_inf_con)*exp(-delta_con*t) + p_inf_con;
    p_dot_con(i) = -delta_con*(p_0_con-p_inf_con)*exp(-delta_con*t);

    
    % error
    e_con(:,i) = J_hat((i-1)*m+1:i*m,:)*th_dot(:,i) - rd_dot(:,i) + k_pos*(r_actual(:,i)-r_desired(:,i)) ...
                 +k_int*(r_aint(:,i)-r_dint(:,i));
    e_cont_norm(:,i) = norm(e_con(:,i));
    % trannsformations
    eta_con(:,i) = e_con(:,i)/p_con(i);
     
    diff_g_xi = diag(exp(xi_u(:,i))./((1+exp(xi_u(:,i))).^2));
    %% control effort
    J_tilde = J_hat((i-1)*m+1:i*m,:)*del*diff_g_xi;
    xi_u_dot = pinv(J_tilde)*(-J_hat((i-1)*m+1:i*m,:)*(nu_neg_dot(:,i)+del_dot*g_xi(:,i))...
                              -J_hat_dot*th_dot(:,i) + rd_ddot(:,i)...
                              -k_pos*(J_hat((i-1)*m+1:i*m,:)*th_dot(:,i) - rd_dot(:,i))...
                              -k_int*(r_actual(:,i) - r_desired(:,i)) ...
                              + eta_con(:,i)*p_dot_con(i)- gamma_c*p_con(i)*diag((1-eta_con(:,i).^2)./2)*log((1+eta_con(:,i))./(1-eta_con(:,i)))); 
    if i<=3 
        xi_u(:,i+1) = xi_u(:,i) + xi_u_dot*dt;
        g_xi(:,i+1) = exp(xi_u(:,i+1))./(1+exp(xi_u(:,i+1))); 
        th_dot(:,i+1) = nu_neg(:,i)+ del*g_xi(:,i+1);
        th(:,i+1) = th(:,i) + th_dot(:,i)*dt;
    else %if i <= 6
         xi_u(:,i+1) = (8/5)*xi_u_dot*dt + (3/5)*xi_u(:,i) + (1/5)*xi_u(:,i-1)+ (1/5)*xi_u(:,i-2);
         g_xi(:,i+1) = exp(xi_u(:,i+1))./(1+exp(xi_u(:,i+1))); 
         th_dot(:,i+1) = nu_neg(:,i) + del*g_xi(:,i+1);
         th(:,i+1) = (8/5)*th_dot(:,i)*dt + (3/5)*th(:,i) + (1/5)*th(:,i-1)+ (1/5)*th(:,i-2);
    
%     else
%         xi_u(:,i+1) = (216/83)*xi_u_dot*dt - (45/83)*xi_u(:,i) + (84/83)*xi_u(:,i-1)+ (82/83)*xi_u(:,i-2) - (27/83)*xi_u(:,i-3) - (21/83)*xi_u(:,i-4) + (10/83)*xi_u(:,i-5);
%         g_xi(:,i+1) = exp(xi_u(:,i+1))./(1+exp(xi_u(:,i+1))); 
%         th_dot(:,i+1) = th_dot_neg + del*g_xi(:,i+1);
%         th(:,i+1) = (216/83)*th_dot(:,i)*dt - (45/83)*th(:,i) + (84/83)*th(:,i-1)+ (82/83)*th(:,i-2) - (27/83)*th(:,i-3) - (21/83)*th(:,i-4) + (10/83)*th(:,i-5);        
     end
    
   th_ddot = del*diff_g_xi*xi_u_dot; 

                     
   error_traj(:,i) = r_actual(:,i) - r_desired(:,i);
   
   %% Observer
   ra_dot(:,i) = J((i-1)*m+1:i*m,:)*(th_dot(:,i)) ;
   ra_ddot(:,i) = J_dot((i-1)*m+1:i*m,:)*(th_dot(:,i)) + J((i-1)*m+1:i*m,:)*del*diff_g_xi*xi_u_dot;
   
   th_dot_inv = pinv((th_dot(:,i)));
%    if all(th_dot(:,i) == 0) 
%        th_dot_inv = zeros(1,n);
%    else
%        th_dot_inv = th_dot(:,i)'/(th_dot(:,i)'*th_dot(:,i));
%    end
    e_obs(:,i) = J_hat((i-1)*m+1:i*m,:)*th_dot(:,i)-ra_dot(:,i);
    e_obs_norm(:,i) = norm(e_obs(:,i));
    e_obs_act = (e_obs(:,i));
%   e_obs_act = bar_lyap(e_obs(:,i),t,norm(e_obs(:,1))+4);
    p_obs(i) = (p_0_obs-p_inf_obs)*exp(-delta_obs*t) + p_inf_obs;
    p_dot_obs(i) = -delta_obs*(p_0_obs-p_inf_obs)*exp(-delta_obs*t);
    eta_obs(:,i) = e_obs(:,i)/p_obs(i);
    % Observer ZNN
    J_hat_dot = (-J_hat((i-1)*m+1:i*m,:)*del*diff_g_xi*xi_u_dot + ra_ddot(:,i)+ eta_obs(:,i)*p_dot_obs(i) - gamma_o*p_obs(i)*diag((1-eta_obs(:,i).^2)./2)*log((1+eta_obs(:,i))./(1-eta_obs(:,i))))*th_dot_inv;
    
    % J_hat for next time step
%     
     if i<=3 
         J_hat((i)*m+1:(i+1)*m,:) = J_hat((i-1)*m+1:i*m,:) + J_hat_dot*dt; %chi update
     else %if i <= 6
          J_hat((i)*m+1:(i+1)*m,:) = (8/5)*J_hat_dot*dt + (3/5)*J_hat((i-1)*m+1:(i)*m,:) + (1/5)*J_hat((i-2)*m+1:(i-1)*m,:)+ (1/5)*J_hat((i-3)*m+1:(i-2)*m,:);
%      else
%          J_hat((i)*m+1:(i+1)*m,:) = (216/83)*J_hat_dot*dt - (45/83)*J_hat((i-1)*m+1:(i)*m,:) + (84/83)*J_hat((i-2)*m+1:(i-1)*m,:) + (82/83)*J_hat((i-3)*m+1:(i-2)*m,:) - (27/83)*J_hat((i-4)*m+1:(i-3)*m,:) - (21/83)*J_hat((i-5)*m+1:(i-4)*m,:) + (10/83)*J_hat((i-6)*m+1:(i-5)*m,:);
    end
   control_eff(:,i) = norm(th_dot(:,i));
   norm_J(i) = norm(J((i-1)*m+1:i*m,:));
   norm_J_dot(i) = norm(J_dot((i-1)*m+1:i*m,:));
   norm_J_hat(i) = norm(J_hat((i-1)*m+1:i*m,:));
   norm_J_hat_dot(i) = norm(J_hat_dot);
   mu(i) = sqrt(det(J_hat((i-1)*m+1:i*m,:)*J_hat((i-1)*m+1:i*m,:)')); 
   r(i) = rank(J_hat((i-1)*m+1:i*m,:));
   pinv_th_dot(:,i) = th_dot_inv';
   norm_pinv_th_dot(:,i) = norm(pinv_th_dot(:,i));
   J_states(i,:) = reshape(J((i-1)*m+1:i*m,:),[1,m*n]);
   J_hat_states(i,:) = reshape(J_hat((i-1)*m+1:i*m,:),[1,m*n]); 
   norm_J_hat_inv(:,i) = norm(pinv(J_hat((i-1)*m+1:i*m,:)));
   norm_J_inv(:,i) = norm(pinv(J((i-1)*m+1:i*m,:)));
   g_xi_norm(:,i) = norm(g_xi(:,i));
   th_ddot_norm(:,i) = norm(th_ddot);
   th_dot_unres(:,i) =  pinv(Jacobian_matrix_kuka(a8,d1,d3,d5,th_unres(:,i)))*rd_dot(:,i);
   th_unres(:,i+1) = th_dot_unres(:,i)*dt + th_unres(:,i);
   error_trag_norm(:,i) = norm(error_traj(:,i));
   norm_jacob_diff(:,i) = norm(J((i-1)*m+1:i*m,:)-J_hat((i-1)*m+1:i*m,:));
end

%% plots
figure(1)
plot(t_span,th_dot(1,1:end-1),t_span,th_dot(2,1:end-1),'c--',t_span,th_dot(3,1:end-1),'y:',t_span,th_dot(4,1:end-1),'g-.',t_span,th_dot(5,1:end-1),'m',t_span,th_dot(6,1:end-1),'k--',t_span,th_dot(7,1:end-1), 'LineWidth', 1.5)
hold on
plot(t_span,th_dot_pos_plot(2,:),'--r',t_span,th_dot_neg_plot(2,:),'--r','LineWidth',1.5)
%title('Variation of the states of system with time')
h2 =legend('$\dot{\theta}_1 $','$\dot{\theta}_2$','$\dot{\theta}_3$','$\dot{\theta}_4$','$\dot{\theta}_5$','$\dot{\theta}_6$','$\dot{\theta}_7$','$\dot{\theta}^+ = 0.3$','$\dot{\theta}^- = -0.3 $',[400 260 0 0]);
set(h2,'Interpreter', 'latex');
xlabel('t(s)')
ylim([-1.05*(max(th_dot_pos)) 1.05*max(th_dot_pos)])
hl = ylabel('$\dot{\theta}_{i}$');
set(hl,'Interpreter', 'latex');

figure(2)
plot(t_span,th(1,1:end-1),t_span,th(2,1:end-1),'c--',t_span,th(3,1:end-1),':',t_span,th(4,1:end-1),'g-.',t_span,th(5,1:end-1),'m',t_span,th(6,1:end-1),'k--',t_span,th(7,1:end-1),'LineWidth',1.5)
hold on
plot(t_span,th_pos_plot(2,:),'--r',t_span,th_neg_plot(2,:),'--r','LineWidth',1.5)
%title('Variation of the states of system with time')
h2 =legend('${\theta}_1 $','${\theta}_2$','${\theta}_3$','${\theta}_4$','${\theta}_5$','${\theta}_6$','${\theta}_7$',[400 260 0 0]);
set(h2,'Interpreter', 'latex');
xlabel('t(s)')
hl = ylabel('${\theta}_{i}$');
ylim([-1.05*(max(th_pos)) 1.05*(max(th_pos))])
set(hl,'Interpreter', 'latex');

figure(3)
plot(t_span(1:(t_end/dt)+1),e_con(1,1:(t_end/dt)+1),t_span(1:(t_end/dt)+1),e_con(2,1:(t_end/dt)+1),'--',t_span(1:(t_end/dt)+1),e_con(3,1:(t_end/dt)+1),'-.')
hold on
title('Variation of controller with time')
plot(t_span(1:(t_end/dt)+1),p_con(1:(t_end/dt)+1),'r--',t_span(1:(t_end/dt)+1),-p_con(1:(t_end/dt)+1),'r--')
hold on
xlabel('t(s)')
ylabel('controller  (m/s)')
legend('e_x','{e_y}','{e_z}','p_u','p_l')
ylabel('controller  (m/s)')
legend('e_x','{e_y}','p_u','p_l')

figure(4)
plot(t_span(1:(t_end/dt)+1),e_obs(1,1:(t_end/dt)+1),t_span(1:(t_end/dt)+1),e_obs(2,1:(t_end/dt)+1),'--',t_span(1:(t_end/dt)+1),e_obs(3,1:(t_end/dt)+1),'-.')
hold on
title('Variation of obs with time')
plot(t_span(1:(t_end/dt)+1),p_obs(1:(t_end/dt)+1),'r--',t_span(1:(t_end/dt)+1),-p_obs(1:(t_end/dt)+1),'r--')
hold on
xlabel('t(s)')
ylabel('observer error (e)  (m/s)')
legend('e_x','{e_y}','{e_z}','p_u','p_l')
legend('e_x','{e_y}','p_u','p_l')



figure(5)
plot(t_span(1:(t_end/dt)+1),error_traj(1,1:(t_end/dt)+1),t_span(1:(t_end/dt)+1),error_traj(2,1:(t_end/dt)+1),'--',t_span(1:(t_end/dt)+1), error_traj(3,1:(t_end/dt)+1),'-.')
hold on
title('Variation of position error with time')
xlabel('t(s)')
ylabel('Position tracking error(\epsilon) (m)')
legend('\epsilon_x','\epsilon_y','\epsilon_z')


figure(6)
plot(t_span, e_cont_norm, t_span, e_obs_norm)
ylabel('||e|||_2')
legend('||e_{cont}|||_2','||e_{obs}|||_2')

figure(7)
% simple tragectory
% simple tragectory
plot3(r_desired(1,:),r_desired(2,:),r_desired(3,:),'b','linewidth',1.5)
hold on

plot3(r_actual(1,:),r_actual(2,:),r_actual(3,:),'r--','linewidth',1.5)
hold on
legend('desired','Traced')


figure(8)
plot(t_span,control_eff)
cont_eff_total = norm(control_eff);
xlabel('time')
ylabel('control effort')
hold on


figure(9)
subplot(1,2,1)
title('Jacobian Norm')
plot(t_span,norm_J, t_span, norm_J_hat)
xlabel('time')
ylabel('Jacobian Norms')
legend('||J||','||J_hat||')
subplot(1,2,2)
title('Difference in Jacobian Norm')
plot(t_span,norm_jacob_diff)
xlabel('time')
ylabel('||J - J_hat||')



figure(11)
plot(t_span, pinv_th_dot)
xlabel('Time(s)')
ylabel('th dot inv')
ylim([-1000 1000])
hold on

figure(12)
plot(t_span, mu)
ylabel('Manipulability Index')
xlabel('Time')
hold on

figure(13)
plot(t_span, r)
ylabel('Rank')
xlabel('Time')
hold on

figure(14)
title('Inverse Norm comparison')
plot(t_span,norm_J_inv, t_span, norm_J_hat_inv)
xlabel('time')
ylabel('Inv Norms')
legend('J inv (Actual)','J hat inv (Estimated)')

figure(15)
title('What Spikes First')
subplot(2,3,1)
plot(t_span, norm_J_hat)
xlabel('Time')
ylabel('||theta dot||')

subplot(2,3,2)
plot(t_span, norm_J_hat_dot)
xlabel('Time')
ylabel('||J Hat dot||')

subplot(2,3,3)
plot(t_span, norm_J_hat_inv)
xlabel('Time')
ylabel('||J hat||')

subplot(2,3,4)
plot(t_span, th_ddot_norm)
xlabel('Time')
ylabel('||theta ddot||')

subplot(2,3,5)
plot(t_span, g_xi_norm)
xlabel('Time')
ylabel('||g xi||')

subplot(2,3,6)
plot(t_span(100:end),norm_pinv_th_dot(100:end))
xlabel('Time')
ylabel('||theta dot inv||')

figure(16)
plot(t_span,th_dot_unres(1,1:end),t_span,th_dot_unres(2,1:end),'c--',t_span,th_dot_unres(3,1:end),'y:',t_span,th_dot_unres(4,1:end),'g-.',t_span,th_dot_unres(5,1:end),t_span,th_dot_unres(6,1:end),t_span,th_dot_unres(7,1:end),'m',t_span,th_dot_pos_plot(1,:),'--r',t_span,th_dot_neg_plot(1,:),'--r','LineWidth',1.5)
%title('Variation of the states of system with time')
h2 =legend('$\dot{\theta}_1 $','$\dot{\theta}_2$','$\dot{\theta}_3$','$\dot{\theta}_4$','$\dot{\theta}_5$','$\dot{\theta}^+ = 0.3$','$\dot{\theta}^- = -0.3 $',[400 260 0 0]);
set(h2,'Interpreter', 'latex');
xlabel('t(s)')

figure(17)
subplot(2,2,1)
plot(t_span,nu_pos)
xlabel('Time')
ylabel('\nu^+')

subplot(2,2,2)
plot(t_span,nu_neg)
xlabel('Time')
ylabel('\nu^-')

subplot(2,2,3)
plot(t_span,nu_pos_dot)
xlabel('Time')
ylabel('\nu dot^+')


subplot(2,2,4)
plot(t_span,nu_neg_dot)
xlabel('Time')
ylabel('\nu dot^-')

figure(18)
subplot(3,2,1)
plot(t_span,th(1,1:end-1), t_span, th_pos_plot(1,:),'r--',t_span, th_neg_plot(1,:),'r--')
legend('\theta_1','\theta_1^+','\theta_1^-')
xlabel('t')
ylabel('\theta_1')
ylim([1.05*th_neg(1) 1.05*th_pos(1)])

subplot(3,2,2)
plot(t_span,th(2,1:end-1), t_span, th_pos_plot(2,:),'r--',t_span, th_neg_plot(2,:),'r--')
legend('\theta_2','\theta_2^+','\theta_2^-')
xlabel('t')
ylabel('\theta_2')
ylim([1.05*th_neg(2) 1.05*th_pos(2)])

subplot(3,2,3)
plot(t_span,th(3,1:end-1), t_span, th_pos_plot(3,:),'r--',t_span, th_neg_plot(3,:),'r--')
legend('\theta_3','\theta_3^+','\theta_3^-')
xlabel('t')
ylabel('\theta_3')
ylim([1.05*th_neg(3) 1.05*th_pos(3)])

subplot(3,2,4)
plot(t_span,th(4,1:end-1), t_span, th_pos_plot(4,:),'r--',t_span, th_neg_plot(4,:),'r--')
legend('\theta_4','\theta_4^+','\theta_4^-')
xlabel('t')
ylabel('\theta_4')
ylim([1.05*th_neg(4) 1.05*th_pos(4)])

subplot(3,2,5)
plot(t_span,th(5,1:end-1), t_span, th_pos_plot(5,:),'r--',t_span, th_neg_plot(5,:),'r--')
legend('\theta_5','\theta_5^+','\theta_5^-')
xlabel('t')
ylabel('\theta_5')
ylim([1.05*th_neg(5) 1.05*th_pos(5)])

subplot(3,2,6)
plot(t_span,th(6,1:end-1), t_span, th_pos_plot(6,:),'r--',t_span, th_neg_plot(6,:),'r--')
legend('\theta_6','\theta_6^+','\theta_6^-')
xlabel('t')
ylabel('\theta_6')
ylim([1.05*th_neg(6) 1.05*th_pos(6)])

figure(19)
subplot(3,2,1)
plot(t_span,th_dot(1,1:end-1), t_span, nu_pos(1,:),'r--',t_span, nu_neg(1,:),'r--')
legend('\theta dot_1','\nu_1^+','\nu_1^-')
xlabel('t')
ylabel('\theta dot_1')
ylim([1.05*nu_neg(1,1) 1.05*nu_pos(1,1)])

subplot(3,2,2)
plot(t_span,th_dot(2,1:end-1), t_span, nu_pos(2,:),'r--',t_span, nu_neg(2,:),'r--')
legend('\theta dot_2','\nu_2^+','\nu_2^-')
xlabel('t')
ylabel('\theta dot_2')
ylim([1.05*nu_neg(2,1) 1.05*nu_pos(2,1)])

subplot(3,2,3)
plot(t_span,th_dot(3,1:end-1), t_span, nu_pos(3,:),'r--',t_span, nu_neg(3,:),'r--')
legend('\theta dot_3','\nu_3^+','\nu_3^-')
xlabel('t')
ylabel('\theta dot_3')
ylim([1.05*nu_neg(3,1) 1.05*nu_pos(3,1)])


subplot(3,2,4)
plot(t_span,th_dot(4,1:end-1), t_span, nu_pos(4,:),'r--',t_span, nu_neg(4,:),'r--')
legend('\theta dot_4','\nu_4^+','\nu_4^-')
xlabel('t')
ylabel('\theta dot_4')
ylim([1.05*nu_neg(4,1) 1.05*nu_pos(4,1)])


subplot(3,2,5)
plot(t_span,th_dot(5,1:end-1), t_span, nu_pos(5,:),'r--',t_span, nu_neg(5,:),'r--')
legend('\theta dot_5','\nu_5^+','\nu_5^-')
xlabel('t')
ylabel('\theta dot_5')
ylim([1.05*nu_neg(5,1) 1.05*nu_pos(5,1)])


subplot(3,2,6)
plot(t_span,th_dot(6,1:end-1), t_span, nu_pos(6,:),'r--',t_span, nu_neg(6,:),'r--')
legend('\theta dot_6','\nu_6^+','\nu_6^-')
xlabel('t')
ylabel('\theta dot_6')
ylim([1.05*nu_neg(6,1) 1.05*nu_pos(6,1)])

figure(20)
semilogy(t_span, error_trag_norm)
ylabel('||\epsilon||')
xlabel('t(s)')
hold on



max_th_1 = max(th(1,:))
min_th_1 = min(th(1,:))


max_th_2 = max(th(2,:))
min_th_2 = min(th(2,:))

max_th_3 = max(th(3,:))
min_th_3 = min(th(3,:))

max_th_4 = max(th(4,:))
min_th_4 = min(th(4,:))

max_th_5 = max(th(5,:))
min_th_5 = min(th(5,:))

max_th_6 = max(th(6,:))
min_th_6 = min(th(6,:))

rel_init_jac_diff = norm_jacob_diff(:,1)/norm_J(:,1)*100
mean_trac_err = norm(error_trag_norm)*sqrt(dt/(m*t_end))
normal_control_eff_new = norm(control_eff)*sqrt(dt/(n*t_end))