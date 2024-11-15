function X = forward_map_5R(L, theta)
    th1 = theta(1);
    th2 = theta(2);
    th3 = theta(3);
    th4 = theta(4);
    th5 = theta(5);
    
    X =L*[cos(th1) + cos(th1+th2) + cos(th1+th2+th3) + cos(th1+th2+th3+th4) + cos(th1+th2+th3+th4+th5);
          sin(th1) + sin(th1+th2) + sin(th1+th2+th3) + sin(th1+th2+th3+th4) + sin(th1+th2+th3+th4+th5)];
 
end