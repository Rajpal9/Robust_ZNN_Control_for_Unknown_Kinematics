%%%%% 2D path for 5R
function [rd, rd_dot, rd_ddot] = path_2D(c,rd_0,t, shape)
% Input arguements
% c - size parameter
% rd_0 -  starting point for the path
% t - time instant
%
    if strcmpi(shape,'adobe') == 1
        rd = [2*c*cos(0.4*pi*t) + 3*c*cos(0.2*pi*t) + rd_0(1) - 5*c;
               2*c*sin(0.4*pi*t) - 3*c*sin(0.2*pi*t) + rd_0(2)];
        rd_dot = [-0.8*pi*c*sin(0.4*pi*t) - 0.6*pi*c*sin(0.2*pi*t);
                   0.8*pi*c*cos(0.4*pi*t) - 0.6*pi*c*cos(0.2*pi*t)];
                   
        rd_ddot = [-0.32*pi^2*c*cos(0.4*pi*t) - 0.12*pi^2*c*cos(0.2*pi*t);
                  -0.32*pi^2*c*sin(0.4*pi*t) + 0.12*pi^2*c*sin(0.2*pi*t)];
              
    elseif strcmpi(shape,'cool_2') == 1
       rd = [2*c*cos(0.2*pi*t) - c*cos(0.8*pi*t) + rd_0(1) - c;
             2*c*sin(0.2*pi*t) - c*sin(0.8*pi*t) + + rd_0(2)];
                   
        rd_dot = [-0.4*c*pi*sin(0.2*pi*t)+ 0.8*c*pi*sin(0.8*pi*t);
                   0.4*c*pi*cos(0.2*pi*t)- 0.8*c*pi*cos(0.8*pi*t)];
                   
        rd_ddot = [-0.08*c*pi^2*cos(0.2*pi*t)+ 0.64*c*pi^2*cos(0.8*pi*t);
                 -0.08*c*pi^2*sin(0.2*pi*t) + 0.64*c*pi^2*sin(0.8*pi*t)];
        
    
    elseif strcmpi(shape,'circle') == 1
       rd =  [c*cos(0.2*pi*t) + rd_0(1)  - c;
              c*sin(0.2*pi*t) + rd_0(2)];
                   
       rd_dot = [-c*0.2*pi*sin(0.2*pi*t);
                 c*0.2*pi*cos(0.2*pi*t)];
                   
       rd_ddot = [-c*((0.2*pi)^2)*cos(0.2*pi*t);
                  -c*((0.2*pi)^2)*sin(0.2*pi*t)];
    
    elseif strcmpi(shape,'cardioid') == 1
       rd = [2*c*cos(0.2*pi*t) - c*cos(0.4*pi*t) + rd_0(1) - c;
             2*c*sin(0.2*pi*t) - c*sin(0.4*pi*t) + rd_0(2)];
                   
       rd_dot = [-2*c*0.2*pi*sin(0.2*pi*t)+c*0.4*pi*sin(0.4*pi*t);
                  2*c*0.2*pi*cos(0.2*pi*t) - c*0.4*pi*cos(0.4*pi*t)];
                   
        rd_ddot = [-2*c*((0.2*pi)^2)*cos(0.2*pi*t)+c*((0.4*pi)^2)*cos(0.4*pi*t);
                   -2*c*((0.2*pi)^2)*sin(0.2*pi*t) + c*((0.4*pi)^2)*sin(0.4*pi*t)];
             
             
    
    elseif strcmpi(shape,'tricuspid') == 1
       rd = [c*cos(0.4*pi*t) + 2*c*cos(0.2*pi*t) + rd_0(1) - 3*c;
             c*sin(0.4*pi*t) - 2*c*sin(0.2*pi*t) + rd_0(2)];
                   
        rd_dot = [-c*0.4*pi*sin(0.4*pi*t)- 2*c*0.2*pi*sin(0.2*pi*t);
                    c*0.4*pi*cos(0.4*pi*t) - 2*c*0.2*pi*cos(0.2*pi*t)];
                   
        rd_ddot = [-c*((0.4*pi)^2)*cos(0.4*pi*t)-2*c*((0.2*pi)^2)*cos(0.2*pi*t);
                   -c*((0.4*pi)^2)*sin(0.4*pi*t) + 2*c*((0.2*pi)^2)*sin(0.2*pi*t)];
            
    
    elseif strcmpi(shape,'star') == 1
       rd = [c*cos(0.6*pi*t) + 2*c*cos(0.4*pi*t) + rd_0(1) - 3*c;
             c*sin(0.6*pi*t) - 2*c*sin(0.4*pi*t) + rd_0(1)];
                   
        rd_dot = [-c*0.6*pi*sin(0.6*pi*t)- 2*c*0.4*pi*sin(0.4*pi*t);
                   c*0.6*pi*cos(0.6*pi*t) - 2*c*0.4*pi*cos(0.4*pi*t)];
                   
        rd_ddot = [-c*((0.6*pi)^2)*cos(0.6*pi*t)-2*c*((0.4*pi)^2)*cos(0.4*pi*t);
                   -c*((0.6*pi)^2)*sin(0.6*pi*t) + 2*c*((0.4*pi)^2)*sin(0.4*pi*t)];
    end
        
        
end