
function handle_connection_MESE_CS_shared(connection)

disp("handle_connection was called.")

setenv('OMP_NUM_THREADS', '10')


%% options and path
load('Param.mat')
output_folder=Param.output_folder;

addpath(genpath(Param.path_to_spm))
addpath(genpath(Param.path_to_EPGX))


%% get data
next_acquisition = @connection.next;

sRO = connection.header.encoding.encodedSpace.matrixSize.x/2;
sPE = connection.header.encoding.encodedSpace.matrixSize.y;
sSE = connection.header.encoding.encodedSpace.matrixSize.z;

nContrast = connection.header.encoding.encodingLimits.contrast.maximum+1;
nSet = connection.header.encoding.encodingLimits.set.maximum+1;
sPE = connection.header.encoding.encodingLimits.kspace_encoding_step_1.maximum+1;

matrixSize = [connection.header.encoding.reconSpace.matrixSize.x...
    connection.header.encoding.reconSpace.matrixSize.y...
    connection.header.encoding.reconSpace.matrixSize.z
    ];

acquisition = next_acquisition(); % Call input function to produce the next acquisition.

nCh = size(acquisition.data.data,2);

kdata=zeros(sRO,nCh,sPE*sSE,nContrast);

Set_idx = acquisition.data.header.set(1:end) + 1;


row = acquisition.data.header.kspace_encode_step_1(1:end) + 1;
col = acquisition.data.header.kspace_encode_step_2(1:end) + 1;

TE_idx = acquisition.data.header.contrast(1:end) + 1;

% store data in kdata
kdata(:,:,sub2ind([sPE,sSE,nContrast],row,col,TE_idx))=acquisition.data.data(:,:,1:end);




%% prepare parameter for reconstruction
kdata = permute(kdata,[1 3 2 4]);
kdata=reshape(kdata,sRO,sPE,sSE,nCh,1,[]); %% buffering the echo according to bart convetion [RO,E1,E2,CHA,MAP,CON]

struct_MESE.TR = connection.header.sequenceParameters.TR;

struct_MESE.calSize=connection.header.userParameters.userParameterDouble(5).value;

%% Estimate coil sensitivity

dimref=[struct_MESE.calSize struct_MESE.calSize struct_MESE.calSize];
func_sens = sprintf('ecalib -r %d %d %d -m 1 -d 3 -c0.1 -k 6:6:6', dimref(1), dimref(2), dimref(3));
sensitivity_coil_map = bart(func_sens,kdata(:,:,:,:,:,1));
disp('---- Sensitivity map estimated ----')

%% Parallel + CS
disp('-------------------------------------------------');
disp('********** bart(pics...)  **********');


im_ec = zeros(size(kdata,1),size(kdata,2),size(kdata,3),nContrast);

% construct basis
TE=connection.header.sequenceParameters.TE;
nT2=200;
nB1=80;
filename='basisfile';
nComp=12;
compl=1;
loga=1;
gadgetron.custom.utils.write_basis(TE,nT2,filename,nComp,nB1,compl,loga);

% reco data
writecfl('kdata',kdata)
writecfl('sensitivity_coil_map',sensitivity_coil_map)
regulCmd='bart pics -i50  -e -d5 -S -R W:7:0:0.001 -B basisfile kdata sensitivity_coil_map output_allec';
system(regulCmd);
basis_final=readcfl('basisfile');
coeff_maps=readcfl('output_allec');
im_ec=squeeze(bart('fmac -s 64', basis_final, coeff_maps));

delete('output_allec.hdr')
delete('output_allec.cfl')
delete('kdata.hdr')
delete('kdata.cfl')
delete('sensitivity_coil_map.hdr')
delete('sensitivity_coil_map.cfl')
info_ec.basisfile.nComp=nComp;
info_ec.basisfile.compl=compl;
info_ec.basisfile.loga=loga;
info_ec.basisfile.nT2=nT2;
info_ec.basisfile.nB1=nB1;



% save info in struct
info_ec.TE=connection.header.sequenceParameters.TE;
info_ec.TR=connection.header.sequenceParameters.TR;
info_ec.FA=connection.header.sequenceParameters.flipAngle_deg;
info_ec.reco=regulCmd;
info_ec.func_sens=func_sens;
info_ec.comment='';




%% create nifti


% Update Output folder
if (isfield(connection.header.studyInformation,'studyDate')) % for offline conversion dat>h5, this field is not present
    full_output_folder= [output_folder '/' connection.header.studyInformation.studyDate '/' connection.header.measurementInformation.protocolName ];

else
    full_output_folder= [output_folder '/Offline_RR_date_' datestr(datetime('today')) '/' connection.header.measurementInformation.protocolName];
end

disp(['Nifti files will be created and saved in : ' full_output_folder])

% create folder if it doesn't exist
if (~exist(full_output_folder))

    full_output_folder=fullfile(full_output_folder,['Process' sprintf('%.02d',1)])
    mkdir(full_output_folder)

else
    ProcDir=cellstr(spm_select('FPList',full_output_folder,'dir','^Pro'));
    full_output_folder=fullfile(full_output_folder,['Process' sprintf('%.02d',length(ProcDir)+1)])
    mkdir(full_output_folder)
end
disp(['create' full_output_folder])


%orientation
if acquisition.data.header.slice_dir(3,1)>0
    Orientation='Tra';
end

if acquisition.data.header.slice_dir(1,1)>0
    Orientation='Sag';
end

if acquisition.data.header.slice_dir(2,1)>0
    Orientation='Cor';
end



if exist('info_ec')
    save(fullfile(full_output_folder,'info_ec.mat'),'info_ec');
end




im=im_ec;
Name='EchoImage';



for k=1:size(im,4)
   

    % create nifti file
    nifti_file_name= [full_output_folder '/' Name '_' sprintf('%.3d',k) '.nii'];
    res=[connection.header.encoding.reconSpace.fieldOfView_mm.x/matrixSize(1) connection.header.encoding.reconSpace.fieldOfView_mm.y/matrixSize(2) connection.header.encoding.reconSpace.fieldOfView_mm.z/matrixSize(3)]
    gadgetron.custom.utils.create_nifti_volume(abs(im(:,:,:,k)),nifti_file_name,reference_header(acquisition),res,matrixSize,Orientation)
     disp("create_and_save_nifti setup...")
end







%% Prepare image to send

first_echo=permute(abs(im_ec(:,:,:,1)),[4, 1, 2, 3]);
scale_factor=4095./max(first_echo(:));

for echo= 1:nContrast
    img_to_send{echo}=scale_factor*(permute(abs(im_ec(:,:,:,echo)),[4, 1, 2, 3]));
end

for echo=1:length(img_to_send)
    image = gadgetron.types.Image.from_data(single(img_to_send{echo}), reference_header(acquisition));
    image.header.image_type = gadgetron.types.Image.MAGNITUDE;
    image = gadgetron.custom.utils.correct_image_siemens(image,acquisition,connection);
    image.header.contrast = echo-1;

    disp("Sending image to client.");
    connection.send(image);
end
end

%% suppport functions
function reference = reference_header(acquisition)
% We pick the first header from the header arrays - we need it to initialize the image meta data.
reference = structfun(@(arr) arr(:, 1)', acquisition.data.header, 'UniformOutput', false);
end



