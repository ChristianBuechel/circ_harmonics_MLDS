function  out = circ_harm(r,a,b,c)
% estimates the geometrical form given by the circular harmonics
% coefficients a, b and c
% r defines #stimuli

%rho = linspace(pi/(r/2),2*pi-pi/(r/2),r-1)';
rho = linspace(0,(2-2/r)*pi,r)';
%rho = linspace(0,2*pi,r)';

desmtx = [];
for i = 1:size(a,1)
   desmtx = [desmtx cos(i*rho)];   
end
for i = 1:size(b,1)
   desmtx = [desmtx sin(i*rho)];   
end

out = c+ desmtx * [a;b];
%rho = [0; rho];
%out = [1; out];
polar(rho,out);
end

