% A script to load vision in the same directory of this file
vpath = which('LoadVision2');
vpath = vpath(1:end-13);

if ~exist('edu.ucsc.neurobiology.vision.io.NeuronFile', 'class')
	javaaddpath([vpath 'Vision.jar']);
end
%configfile = [vpath 'config_base.xml'];
configfile = [vpath 'config.xml'];

%defaultsfile = [vpath 'cortexth6_angled_allspike.xml'];
defaultsfile = [vpath 'invivo_v4.xml'];

