function [fits ] = fit_all
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
n_faces = 8;
subj    = [1:12 14 16:36]; 
sess    = [2 4];
for s = 1:size(sess,2)
fits{s}    = zeros(n_faces,size(subj,2));
for i=1:size(subj,2)
    disp(['doing ' num2str(subj(i)) ' session ' num2str(sess(s))]);
    [para, coord] = circ_harm_MLDS(subj(i),sess(s));
    fits{s}(:,i)  = coord;
end

end

