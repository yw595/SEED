baseDir = '/home/ubuntu/MATLAB/SEED';
modelsDir = '/home/ubuntu/MATLAB/SEED/input/SEEDModels';
GreenblumDir = '/home/ubuntu/MATLAB/SEED/input/GreenblumData';
FI = fopen([GreenblumDir filesep 'GreenblumECs.txt']);
GreenblumEC = textscan(FI,'%s\n');
fclose(FI);
GreenblumEC = GreenblumEC{1};
if ~exist('CobraLPSolver','var')
    initCobraToolbox;
end
