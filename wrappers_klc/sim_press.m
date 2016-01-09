function [out_posxy] = sim_press(spin_sys, more_params)

% Jamie Near, McGill University 2015.
% 
% USAGE:
% This script is run simply by editing the input parameters and then
% clicking "Run".
% 
% DESCRIPTION:
% This script simulates a PRESS experiment with fully shaped refocusing 
% pulses.  Phase cycling of refocusing pulses is performed.  Furthermore, 
% simulations are run at various locations in space to account for the 
% within-voxel spatial variation of the metabolite signal.  Summation 
% across phase cycles and spatial positions is performed.  As a result of 
% the phase cycling and spatially resolved simulations, this code takes a 
% long time to run.  Therefore, the MATLAB parallel computing toolbox 
% (parfor loop) was used to accelerate the siumulations.  Acceleration 
% is currently performed in the direction of the slice selective pulse along
% the x-direction, but this can be changed.  Up to a factor of 12 acceleration
% can be achieved using this approach.  To enable the use of the MATLAB
% parallel computing toolbox, initialize the multiple worked nodes using
% "matlabpool size X" where "X" is the number of available processing
% nodes.  If the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.
% 
% INPUTS:
% To run this script, edit the parameters below as desired and then click
% "run":
% refocWaveform     = name of refocusing pulse waveform.
% refTp             = duration of refocusing pulses[ms]
% Bfield            = Magnetic field strength in [T]
% Npts              = number of spectral points
% sw                = spectral width [Hz]
% Bfield            = magnetic field strength [Tesla]
% lw                = linewidth of the output spectrum [Hz]
% thkX              = slice thickness of x refocusing pulse [cm]
% thkY              = slice thickness of y refocusing pulse [cm]
% x                 = vector of X positions to simulate [cm]
% y                 = vector of y positions to simulate [cm]
% taus              = vector of pulse sequence timings  [ms]
% spinSys           = spin system to simulate 
% refPhCyc1         = vector of phase cycling steps for 1st refocusing pulse [degrees]
% refPhCyc2         = vector of phase cycling steps for 2nd refocusing pulse [degrees]

% ************INPUT PARAMETERS**********************************
refocWaveform='sampleRefocPulse.pta'; %name of refocusing pulse waveform.
refTp=2.5; %duration of refocusing pulses[ms]
Npts=8192; %number of spectral points
sw=2000; %spectral width [Hz]
lw=8; %linewidth of the output spectrum [Hz]
Bfield=3; %Magnetic field strength in [T]
thkX=1.66; %slice thickness of x refocusing pulse [cm]
thkY=1.66; %slice thickness of y refocusing pulse [cm]
%x=linspace(-1,1,32); %X positions to simulate [cm]
%y=linspace(-1,1,32); %y positions to simulate [cm]
x=0;
y=0; 
TE = 140;
TE1 = 13.4;
pulse_type = 'shaped';
%timing of the pulse sequence [ms]
centreFreq=3.0; %Centre frequency of MR spectrum [ppm]
refPhCyc1=[0,90]; %phase cycling steps for 1st refocusing pulse [degrees]
refPhCyc2=[0,90]; %phase cycling steps for 2ndtau2 refocusing pulse [degrees]

run check_param_fields

TE2=TE-TE1;

if strcmp(pulse_type, 'ideal')
    pulse_seq = @sim_press_ideal;
elseif strcmp(pulse_type, 'shaped')
    pulse_seq = @sim_press_shaped;
end

% ************END OF INPUT PARAMETERS**********************************

%Load RF waveform
refRF=io_loadRFwaveform(refocWaveform,'ref',0);

gamma=42577000; %gyromagnetic ratio

%Load spin systems
load spinSystems
sys=eval(['sys' spin_sys]);

%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
refRF=rf_resample(refRF,100);

Gx=(refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
Gy=(refRF.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]

[DX,DY]=meshgrid(x,y);

%n=1;
%totalIters=length(x)*length(y)*length(editPhCyc1)*length(editPhCyc2)*length(refPhCyc1)*length(refPhCyc2);

%Initialize structures:
out_posxy_rpc=cell(length(x),length(y),length(refPhCyc1),length(refPhCyc2));
out_posxy=cell(length(x),length(y));
out=struct([]);

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

for X=1:length(x);
%parfor X=1:length(x);
    for Y=1:length(y);
        for RP1=1:length(refPhCyc1)
            for RP2=1:length(refPhCyc2)
                disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) ', '...
                    'Y-position ' num2str(Y) ' of ' num2str(length(y)) ', '...
                    'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) ', '...
                    'Second Refoc phase cycle ' num2str(RP2) ' of ' num2str(length(refPhCyc2)) '!!!']);
                out_posxy_rpc{X}{Y}{RP1}{RP2}=pulse_seq(Npts,sw,Bfield,lw,sys,TE1,TE2,...
                    refRF,refTp,x(X),y(Y),Gx,Gy,refPhCyc1(RP1),refPhCyc2(RP2));
                
                if RP1==1 && RP2==1
                    out_posxy{X}{Y}=out_posxy_rpc{X}{Y}{RP1}{RP2};
                else
                    out_posxy{X}{Y}=op_addScans(out_posxy{X}{Y},out_posxy_rpc{X}{Y}{RP1}{RP2},xor(RP1==length(refPhCyc1),RP2==length(refPhCyc2)));
                end
                out_posxy{X}{Y}.x = x;
            end %end of 1st refocusing phase cycle loop
        end %end of 2nd refocusing phase cycle loop.
        
        out=op_addScans(out,out_posxy{X}{Y});
        
    end %end of spatial loop (parfor) in y direction.
end %end of spatial loop (parfor) in x direction.

end




       