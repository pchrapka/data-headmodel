%% create_icbm152

do_checks = false;

%% create mri
%%%%%%%%%%%%%%%%%%%%%%%%%

% load MRI data for ICMB 152
mri = ft_read_mri('MNI152_T1_0.5mm.nii');
mri.coordsys = 'mni';

if do_checks
    % check fiducial locations
    cfg = [];
    cfg.location = [0 84 -43]; %nas
    ft_sourceplot(cfg,mri);
    
    cfg.location = [-75.09 -19.49 -47.98]; %lpa
    ft_sourceplot(cfg,mri);
    
    cfg.location = [76 -19.45 -47.7]; %rpa
    ft_sourceplot(cfg,mri);
end

% add fiducials to mri header
% MNI coordinates are taken from: Cutini S, Scatturin P, Zorzi M (2011): A
% new method based on ICBM152 head surface for probe placement in
% multichannel fNIRS

fiducials = [];
fiducials.nas = [0 84 -43];
fiducials.lpa = [-75.09 -19.49 -47.98];
fiducials.rpa = [76 -19.45 -47.7];
fid_names = fieldnames(fiducials);

% transform from MNI to original anatomical coordinate space
transform = mri.transform;
R_inv = inv(transform(1:3,1:3));
d = transform(1:3,4);
inv_transform = [R_inv -R_inv*d; 0 0 0 1];

fid_anatomical = [];
for i=1:length(fid_names)
    field = fid_names{i};
    fid_anatomical.(field) = ft_warp_apply(inv_transform, fiducials.(field), 'homogenous');
    
    % check
    if do_checks
        % they should be equal
        fid_mni = ft_warp_apply(transform, fid_anatomical.(field), 'homogenous');
        fprintf('Fiducial: %s\nMNI\n',field);
        disp(fid_mni);
        fprintf('Original\n');
        disp(fiducials.(field));
    end
    
end

mri.hdr.fiducials.mri = fid_anatomical;
mri.hdr.fiducials.head = fid_anatomical;

% save mri
save('icbm152_mri.mat','mri');

%% segment volume
%%%%%%%%%%%%%%%%%%%%%%%%%

% Source: http://www.agricolab.de/template-headmodel-for-fieldtrip-eeg-source-reconstruction-based-on-icbm152/

cfg = [];
cfg.brainthreshold = 0.5;
cfg.scalpthreshold = 0.15;
cfg.downsample = 1;
cfg.output = {'brain','skull','scalp'};

seg = ft_volumesegment(cfg, mri);

save('icbm152_seg.mat','seg');


%% prepare mesh
%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = {'brain','skull','scalp'};
cfg.numvertices = [2000, 2000, 1000];

mesh = ft_prepare_mesh(cfg,seg);

save('icbm152_mesh.mat','mesh');
