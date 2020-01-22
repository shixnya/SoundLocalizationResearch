function subasdf = ASDFSubsample(asdf, subIndex)
% subasdf = ASDFsubsample(asdf, subIndex)
%
%    asdf - {nNeu+2,1} ASDF to subsample from
%    subIndex - (subNeu,1) indices of subsampled neurons
%
% Returns:
%    subasdf - {subNeu+2,1} subsampled ASDF
%
% Description :
%
%
% Example :
%    
%
% Author   : Shinya Ito
%            Indiana University
%
% Modified on 2/4/2010
% Modified on 6/26/2012  l.24 'length' is changed to 'nnz'
%                        Suggenstion by Ben Nicholson

subasdf = asdf(subIndex);

subNeu = nnz(subIndex);

subasdf{end+1} = asdf{end-1};
subasdf{end+1} = [subNeu, asdf{end}(2)]; % new [nNeu duration]
