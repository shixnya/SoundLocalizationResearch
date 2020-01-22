function visionObject = LoadVisionFiles(filename)

if isfullpath(filename) % fullpath
    fullpath = filename;
else
    fullpath = [pwd filesep filename];
end

[path name ext] = fileparts(fullpath);

if ~isempty(strfind(name, '*')) % there is * in string
    exname = getWildName(fullpath);
    fullpath = [path filesep exname];
end
    

if exist(fullpath, 'file');
    switch ext
        case '.spikes'
            visionObject = edu.ucsc.neurobiology.vision.io.SpikeFile(fullpath);
        case '.neurons'
            visionObject = edu.ucsc.neurobiology.vision.io.NeuronFile(fullpath);
        case '.neurons-raw'
            visionObject = edu.ucsc.neurobiology.vision.io.NeuronFile(fullpath);
        case '.params'
            visionObject = edu.ucsc.neurobiology.vision.io.ParametersFile(fullpath);
        case '.bin'
            visionObject = edu.ucsc.neurobiology.vision.io.RawDataFile(fullpath);
        case '.prj'
            visionObject = edu.ucsc.neurobiology.vision.io.ProjectionsFile(fullpath);
        case '.ei'
            visionObject = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(fullpath);
        case '.model'
            visionObject = edu.ucsc.neurobiology.vision.io.ClusteringModelFile(fullpath);
        case '.sta'
            visionObject = edu.ucsc.neurobiology.vision.io.STAFile(fullpath);
        otherwise % assuming raw data
            try
                visionObject = edu.ucsc.neurobiology.vision.io.RawDataFile(fullpath);
            catch
                fullpath
                warning(['LoadVisionFiles: Unknown vision file type ' ext ' is specified. returning 0...']);
                visionObject = 0;
            end
    end
else
    error(['LoadVisionFiles: File does not exist. Check the path. Given Path: ' fullpath]);
end






        

