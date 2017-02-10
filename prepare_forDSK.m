%% prepare_forDSK 
% a Matlab routine that prepare DEM file over a specified  area to be 
%  ingested into DSK toolkit provided by the NAIF Spice toolkit.
%
%
%  /!\ The routine calls dem2dsk and make_dskkernel routines 
%
%
% A. Lucas (dralucas@geophysx.edu.eu.org), 2017
%% Clearing matlab env.
clc;close all hidden;clear all
%% Gloabal variables declaration
global name dispfigure savefigure dlkernel runmkdsk
%% USER INPUTS
filein='ldem_16.lbl';        % DEM, should be in PDS-3 format (for now)
ppd='16';                    % Pixel per degree sampling, it's used for naming convention of the output files.
name={'Hermite'};            % Location name, can ba liste
location = {[86.17 266.6]};  % Location coordinates in decimal degrees [lat,lon] which must fit the DEM file coordinates system,
search_radius={[2.5 25]};    % Search in decimal degrees [Delta_latitude Delta_longitude]; so adjustment can be made for high latitude areas.
dispfigure=false;            % Save figure "True/False"
savefigure=false;            % Save figure "True/False", could take a while
dlkernel=false;              % Download the kernels "True/False", could take a while
runmkdsk=true;
%% Loop over location areas provided - DO NOT EDIT BELOW !!!
for kk=1:numel(name)
    DEG = (search_radius{kk}) ;
     
    latrange = [location{kk}(1)-DEG(1) location{kk}(1)+DEG(1)];
    latrange(latrange<-90) = 90;
    latrange(latrange>90) = 90;
    
    if latrange(1)>latrange(2)
        tmp = latrange(2);
        latrange(2) = latrange(1);
        latrange(1) = tmp;
        clear tmp
    end
    
    if latrange(1)==latrange(2)
        disp('!!! Error: no range in latitude.')
        return;
    end
       
    lonrange = [location{kk}(2)-DEG(2) location{kk}(2)+DEG(2)];   
    
    if lonrange(1)>lonrange(2)
        tmp = lonrange(2);
        lonrange(2) = lonrange(1);
        lonrange(1) = tmp;
        clear tmp
    end
    
    if lonrange(1)==lonrange(2)
        disp('!!! Error: no range in longitude.')
        return;
    end
    
    fileout=[name{kk},'_',ppd,'ppd_plate.tab'];
    [tri,dem,lbl]=dem2dsk(filein,fileout,latrange,lonrange);
    clear tri dem lbl
end
clear all;