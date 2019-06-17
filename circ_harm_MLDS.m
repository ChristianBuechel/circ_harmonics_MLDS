function [fit plotme ll] = circ_harm_MLDS(sub, sess, order)
% use circular harmonics to fit Selim's data
% Detailed explanation goes here

n_faces   = 8;
nb_cos    = order; %# basis functions cosine
nb_sin    = order; %# basis functions sine
% construct all angles
rho       = linspace(0,(2-2/n_faces)*pi,n_faces)';

desmtx = [];
for i = 1:nb_cos
   desmtx = [desmtx cos(i*rho(2:end))];   
   %desmtx = [desmtx cos(i*rho)];   
end
for i = 1:nb_sin
   desmtx = [desmtx sin(i*rho(2:end))];   
   %desmtx = [desmtx sin(i*rho)];   
end

base      = 'c:\Users\buechel\Data\lukas\testData\';
data_file = [base sprintf('sub%3.3d',sub) filesep 'phase' sprintf('%2.2d',sess) filesep 'quadruplet\data.mat'];

global data;
data = load(data_file);
data.out.angles    = rho(data.out.Sequence);
data.out.responses = ~data.out.Response;
data.out.desmtx    = desmtx;
data.out.nb_cos    = nb_cos;
data.out.nb_sin    = nb_sin;
data.out.rho       = rho;


x_init = [zeros(1,nb_cos) zeros(1,nb_sin) 1 1]; %cosines sines constant sd(Ncdf)
%x_init = [randn(1,nb_cos) randn(1,nb_sin) 3+randn(1) 10+randn(1)]; %cosines sines constant sd(Ncdf)
options.MaxFunEvals = 10000;
fit = fminsearch(@negll_faces,x_init,options,0); %optimize

[ll, plotme] = negll_faces(fit,1);

function [ll allcoord] = negll_faces(x,plot)
global data;
% sort all parameters
b_cos  = x(1:data.out.nb_cos)';
b_sin  = x(data.out.nb_cos+1:data.out.nb_cos+data.out.nb_sin)';
b_cons = x(data.out.nb_cos+data.out.nb_sin+1);
b_sd   = x(data.out.nb_cos+data.out.nb_sin+2);

m_radii = b_cons + data.out.desmtx * [b_cos;b_sin]; %estimate radii based on spherical harmonics
m_radii = [1;m_radii]; % radius for face = 1
%m_radii = m_radii - mean(m_radii)+1;
radii   = m_radii(data.out.Sequence);

% now calculate distance of each pair
d_r(:,1)   = diff(radii(:,[1 2]),1,2);
d_r(:,2)   = diff(radii(:,[3 4]),1,2);
d_rho(:,1) = diff(data.out.angles(:,[1 2]),1,2);
d_rho(:,2) = diff(data.out.angles(:,[3 4]),1,2);
dist(:,1)  = sqrt(sum(radii(:,1:2).^2,2)-2.*prod(radii(:,1:2),2).*cos(d_rho(:,1)));
dist(:,2)  = sqrt(sum(radii(:,3:4).^2,2)-2.*prod(radii(:,3:4),2).*cos(d_rho(:,2)));
d_dist     = diff(dist,1,2);

p = spm_Ncdf(d_dist,0,b_sd).^(data.out.responses) .* (1-spm_Ncdf(d_dist,0,b_sd)).^(1-data.out.responses);
ll = -sum(log(p));
if plot
    allcoord = m_radii; 
    %clf;
%     subplot(2,1,1);
%     xx=linspace(-1.5,1.5,50);
%     plot(xx,spm_Ncdf(xx,0,b_sd),'b-')
%     hold on
%     plot(d_dist,data.out.responses,'rx')
%     subplot(2,1,2);
    polar(data.out.rho,m_radii,'r-');
    view([90 -90]);drawnow;
end
