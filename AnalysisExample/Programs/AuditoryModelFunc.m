function dist = AuditoryModelFunc(p, x)
% A function that defines auditory responses of the neurons
% get 
% x{3} is ILD_mini.

w1 = p(1); % weight for spectral cue
w2 = p(2); % weight for ILD
bl = p(3); % baseline firing rate
wid = p(4); % width of the spectral cue Kent (kappa)
%spp = p(5); % parameter for softplus

pk = [wid, 0, 0, 0, 0, w1, bl];

dist = max(kentdist(pk, x) + w2 * x{3}, 0);
%dist = softplus(kentdist(pk, x) + w2 * x{3}, spp);

