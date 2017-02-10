%% make_dskkernels
% a Matlab subroutine that prepare kernel lists that is requiered for
%  the DSK toolkit provided by the NAIF Spice toolkit.
%
%
%  /!\ The routine is constantly evolving, no stable version available
%
%
% A. Lucas (dralucas@geophysx.edu.eu.org), 2017

function [kernels]=make_dskkernels(body)
global dlkernel
switch body
    case 'moon'
        kernels.list ={'moon_pa_de421_1900-2050.bpc','moon_080317.tf',...
            'moon_assoc_me.tf','pck00010.tpc','naif0012.tls','de430.bsp'};
        kernels.SURFACE_NAME='MOON';
        kernels.REF_FRAME_NAME='IAU_MOON';
end

if dlkernel
    %dl the required kernels
    get_kernels(kernels);
end

end




function get_kernels(kernels)

url='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/';
disp(' ')
disp('[GET SPICE KERNELS]')
for kk=1:numel(kernels.list)
    [p,n,e] = fileparts(kernels.list{kk});
    switch e
        case '.tf'
            urlfile=[url,'fk/satellites/',kernels.list{kk}];
        case '.tls'
            urlfile=[url,'lsk/',kernels.list{kk}];
        case '.bpc'
            urlfile=[url,'pck/',kernels.list{kk}];
        case '.tpc'
            urlfile=[url,'pck/',kernels.list{kk}];
        case '.bsp'
            urlfile=[url,'spk/planets/',kernels.list{kk}];
    end
    
    disp(['  *Get: ',kernels.list{kk}]);
    %urlwrite(urlfile,kernels.list{kk});
    if isunix
        if ismac
          unix(['/usr/local/bin/wget -c  --no-check-certificate ',urlfile]) % We assume that wget is installed from Homebrew/Macports...
        else 
           unix(['wget -c  --no-check-certificate ',urlfile])
        end
    else
        disp('!!! Error: Only Un*X/POSIX systems are supported.')
    end
    
end
end
