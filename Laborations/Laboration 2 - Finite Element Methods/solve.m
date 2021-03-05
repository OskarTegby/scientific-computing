%=======================================================================
% Numerical Solver
%=======================================================================

function [solution, tvec]=solve(dx,dt,method,nframes,Tfinal, ...
                                u1_left,u2_left,u3_left,u1_right, ...
                                u2_right,u3_right);  
disp('COMPUTING...');
if (method==1)
  [solution, tvec]=leapfrog(dx,dt,nframes,Tfinal,u1_left,u2_left, ...
                            u3_left,u1_right,u2_right,u3_right); 
else
  [solution, tvec]=rusanov(dx,dt,nframes,Tfinal,u1_left,u2_left, ...
                           u3_left,u1_right,u2_right,u3_right); 
end
disp('READY');
