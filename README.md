# Example code to reconstruct and retrospectively correct for motion. 

To be used with data acquired with the 3D cartesian Multi-Echo Spin echo sequence.

*Nadege Corbin, 
Centre de Resonance Magnetique des Systemes Biologiques 
UMR 5536 CNRS/Bordeaux University 
2024.04.16*


# Requirements


This reconstruction code requires multiple softwares and toolbox

- MATLAB (tested with R2019b) 
- Gadgetron (https://github.com/gadgetron/gadgetron)
- siemens_to_ismrmrd (https://github.com/ismrmrd/siemens_to_ismrmrd)
- SPM (https://github.com/spm/)
- EPG-X (https://github.com/mriphysics/EPG-X)

# Step 1: copy files in the right place

- create the folder `[path to gadgetron]/+gadgetron/+custom/+MESE_CS`
- copy `handle_connection_MESE_CS_shared.m` there
- if it does not exist, create the folder ```
[path to gadgetron]/+gadgetron/+custom/+utils ```
- copy `write_basis_shared.m` there 

- copy `MESE_CS_shared.xml  into `/usr/local/share/gadgetron/config`


# Step 2: Update path 

update the path to SPM, EPG-X and sdc3_nrz_11aug  in the matlab structure `Param.mat`
update the poutput_folder in in the matlab structure `Param.mat`

# Step 3 : Convert raw data (.dat) in Ismrmrd format (.h5)

siemens_to_ismrmrd -f [dat_file_name] -o [h5_file_name] -x IsmrmrdParameterMap_Siemens.xsl -z2 


# Step 4: start gadgetron
Open a terminal 
go to the "Reco_MESE_CS" folder
start gadgetron with the  command: `gadgetron -p 9002`

# Step 5: Reconstruct Data 
Open a second terminal:
go to the "Reco_MESE_CS" folder
reconstruct the data with the following command
`gadgetron_ismrmrd_client -p 9002 -f [h5_file_name] -c MESE_CS_shared.xml `





