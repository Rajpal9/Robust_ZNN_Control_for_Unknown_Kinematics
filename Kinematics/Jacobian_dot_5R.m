function J_dot = Jacobian_dot_5R(L,theta,theta_dot)
    th1 = theta(1);
    th2 = theta(2);
    th3 = theta(3);
    th4 = theta(4);
    th5 = theta(5);
    
    th1_dot = theta_dot(1);
    th2_dot = theta_dot(2);
    th3_dot = theta_dot(3);
    th4_dot = theta_dot(4);
    th5_dot = theta_dot(5);
    
     J_dot = L*[-(cos(th1)*(th1_dot) + cos(th1+th2)*(th1_dot+th2_dot) + cos(th1+th2+th3)*(th1_dot+th2_dot+th3_dot)+cos(th1+th2+th3+th4)*(th1_dot+th2_dot+th3_dot+th4_dot)+cos(th1+th2+th3+th4+th5)*(th1_dot+th2_dot+th3_dot+th4_dot+th5_dot)), -(cos(th1+th2)*(th1_dot+th2_dot) + cos(th1+th2+th3)*(th1_dot+th2_dot+th3_dot)+cos(th1+th2+th3+th4)*(th1_dot+th2_dot+th3_dot+th4_dot)+cos(th1+th2+th3+th4+th5)*(th1_dot+th2_dot+th3_dot+th4_dot+th5_dot)), -(cos(th1+th2+th3)*(th1_dot+th2_dot+th3_dot)+cos(th1+th2+th3+th4)*(th1_dot+th2_dot+th3_dot+th4_dot)+cos(th1+th2+th3+th4+th5)*(th1_dot+th2_dot+th3_dot+th4_dot+th5_dot)),-(cos(th1+th2+th3+th4)*(th1_dot+th2_dot+th3_dot+th4_dot)+cos(th1+th2+th3+th4+th5)*(th1_dot+th2_dot+th3_dot+th4_dot+th5_dot)),-(cos(th1+th2+th3+th4+th5)*(th1_dot+th2_dot+th3_dot+th4_dot+th5_dot));
             -(sin(th1)*(th1_dot) + sin(th1+th2)*(th1_dot+th2_dot) + sin(th1+th2+th3)*(th1_dot+th2_dot+th3_dot)+sin(th1+th2+th3+th4)*(th1_dot+th2_dot+th3_dot+th4_dot)+sin(th1+th2+th3+th4+th5)*(th1_dot+th2_dot+th3_dot+th4_dot+th5_dot)), -(sin(th1+th2)*(th1_dot+th2_dot) + sin(th1+th2+th3)*(th1_dot+th2_dot+th3_dot)+sin(th1+th2+th3+th4)*(th1_dot+th2_dot+th3_dot+th4_dot)+sin(th1+th2+th3+th4+th5)*(th1_dot+th2_dot+th3_dot+th4_dot+th5_dot)), -(sin(th1+th2+th3)*(th1_dot+th2_dot+th3_dot)+sin(th1+th2+th3+th4)*(th1_dot+th2_dot+th3_dot+th4_dot)+sin(th1+th2+th3+th4+th5)*(th1_dot+th2_dot+th3_dot+th4_dot+th5_dot)),-(sin(th1+th2+th3+th4)*(th1_dot+th2_dot+th3_dot+th4_dot)+sin(th1+th2+th3+th4+th5)*(th1_dot+th2_dot+th3_dot+th4_dot+th5_dot)),-(sin(th1+th2+th3+th4+th5)*(th1_dot+th2_dot+th3_dot+th4_dot+th5_dot))]; 
    end