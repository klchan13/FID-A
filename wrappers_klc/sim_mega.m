function [spect] = sim_mega(spin_sys, more_params)

% run_simMegaPressShaped.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% This script is run simply by editing the input parameters and then
% clicking "Run".
% 
% DESCRIPTION:
% This script simulates a MEGA-PRESS experiment with fully shaped editing 
% and refocusing pulses.  Phase cycling of both the editing and refocusing
% pulses is performed.  Furthermore, simulations are run at various
% locations in space to account for the within-voxel spatial variation of
% the GABA signal.  Summation across phase cycles and spatial positions is
% performed.  As a result of the phase cycling and spatially resolved simulations, 
% this code takes a long time to run.  Therefore, the MATLAB parallel computing
% toolbox (parfor loop) was used to accelerate the siumulations.  Accelration 
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
% editWaveform      = name of editing pulse waveform.
% editOnFreq        = freqeucny of edit on pulse[ppm]
% editOffFreq       = frequency of edit off pulse[ppm]
% refTp             = duration of refocusing pulses[ms]
% editTp            = duration of editing pulses[ms]
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
% editPhCyc1        = vector of phase cycling steps for 1st editing pulse [degrees]
% editPhCyc2        = vector of phase cycling steps for 2nd editing pulse [degrees]
% refPhCyc1         = vector of phase cycling steps for 1st refocusing pulse [degrees]
% refPhCyc2         = vector of phase cycling steps for 2nd refocusing pulse [degrees]

% ************INPUT PARAMETERS**********************************
refocWaveform='sampleRefocPulse.pta'; %name of refocusing pulse waveform.
% refocWaveform1='sampleRefocPulse.pta'; %name of refocusing pulse waveform.
% refocWaveform2='sampleRefocPulse.pta'; %name of refocusing pulse waveform.
editWaveform='sampleEditPulse.pta'; %name of editing pulse waveform.
editOnFreq=4.1; %freqeucny of edit on pulse[ppm]
editOffFreq=10; %frequency of edit off pulse[pp m]
refTp=6.91;%duration of refocusing pulses[ms]
editTp=14;
Npts=8192; %number of spectral points
sw=2000; %spectral width [Hz]
lw=2.5; %linewidth of the output spectrum [Hz]
Bfield=3; %Magnetic field strength in [T]
thkX=3; %slice thickness of x refocusing pulse [cm]
thkY=3; %slice thickness of y refocusing pulse [cm]
nx = 19;
ny = 19;
%x=-1.8:(3.6/(nx-1)):1.8; %X positions to simulate [cm]
%y=-1.8:(3.6/(ny-1)):1.8; %y positions to simulate [cm]
x=0;
y=0;
TE = 144.0;
TE1 = 13.4;
pulse_type = 'shaped';
edit_angle_on = 180;
edit_angle_off = 180;
  
centreFreq=3.0; %Centre frequency of MR spectrum [ppm]
editPhCyc1=[0 90]; %phase cycling steps for 1st editing pulse [degrees]
editPhCyc2=[0 90];% 180 270]; %phase cycling steps for 2nd editing pulse [degrees]
refPhCyc1=[0,90]; %phase cycling steps for 1st refocusing pulse [degrees]
refPhCyc2=[0,90]; %phase cycling steps for 2nd refocusing pulse [degrees]

run check_param_fields

TE2=TE-TE1;%TE1 in ms

taus=[TE1/2 - refTp/2,...
      TE2/4 - refTp/2- editTp/2,...
      TE2/4  + TE1/2 - refTp/2 - editTp/2,...
      TE2/4 - refTp/2 - editTp/2,...
      TE2/4 - editTp/2]; %timing of the pulse sequence [ms]
  
if strcmp(pulse_type, 'edit_ideal')
    pulse_seq = @sim_mega_edit_ideal;
elseif strcmp(pulse_type, 'refoc_ideal')
    pulse_seq = @sim_mega_refoc_ideal;
elseif strcmp(pulse_type, 'shaped')
    pulse_seq = @sim_mega_shaped;
elseif strcmp(pulse_type, 'ideal')
    pulse_seq = @sim_mega_ideal;
end

% ************END OF INPUT PARAMETERS**********************************

%Load RF waveforms
refRF=io_loadRFwaveform(refocWaveform,'ref',0);
%refRF1=io_loadRFwaveform(refocWaveform1,'ref',0);
%refRF2=io_loadRFwaveform(refocWaveform2,'ref',0);
% refRF.tbw = 38.280; %9.32;
% refRF.tw1 = 9.932;
editRF=io_loadRFwaveform(editWaveform,'inv',0);
%editRF.waveform(:,3) = editRF.waveform(:,3)*(0.010/100);
% editRF.tw1 = 0.914;
%editRF.tbw = 1.5;

gamma=42577000; %gyromagnetic ratio

%Load spin systems
load spinSystems
sys=eval(['sys' spin_sys]);

%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
%refRF=rf_resample(refRF,128); % Comment this for FMREF07

%This is the step where the editing pulse waveform (initially a pulse with 
%zero-frequency) is frequency shifted to produce and edit-on and an
%edit-off pulse;

% editRF.waveform = interp1(editRF.waveform,1:0.25:101, 'spline');
% editRF.waveform(:,3) = editRF.waveform(:,3)*(0.010/400);

editRFon=rf_freqshift(editRF,editTp,(centreFreq-editOnFreq)*Bfield*gamma/1e6);
editRFoff=rf_freqshift(editRF,editTp,(centreFreq-editOffFreq)*Bfield*gamma/1e6);

Gx=(refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
Gy=(refRF.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
% Gx1=(refRF1.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
% Gy1=(refRF1.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]

% Gx2=(refRF2.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
% Gy2=(refRF2.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]

[DX,DY]=meshgrid(x,y);

%n=1;
%totalIters=length(x)*length(y)*length(editPhCyc1)*length(editPhCyc2)*length(refPhCyc1)*length(refPhCyc2);

%Initialize structures:
outON_posxy_epc_rpc=cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2),length(refPhCyc1),length(refPhCyc2));
outOFF_posxy_epc_rpc=cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2),length(refPhCyc1),length(refPhCyc2));
outON_posxy_epc=cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2));
outOFF_posxy_epc=cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2));
outON_posxy=cell(length(x),length(y));
outOFF_posxy=cell(length(x),length(y));
outON=struct([]);
outOFF=struct([]);


%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%parfor X=1:length(x);
for X=1:length(x);
    for Y=1:length(y);
        for EP1=1:length(editPhCyc1)
            for EP2=1:length(editPhCyc2)
                for RP1=1:length(refPhCyc1)
                    for RP2=1:length(refPhCyc2)
                            
                        disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) ', '...
                            'Y-position ' num2str(Y) ' of ' num2str(length(y)) ', '...
                            'First Edit phase cycle ' num2str(EP1) ' of ' num2str(length(editPhCyc1)) ','...
                            'Second Edit phase cycle ' num2str(EP2) ' of ' num2str(length(editPhCyc2)) ', '...
                            'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) ', '...
                            'Second Refoc phase cycle ' num2str(RP2) ' of ' num2str(length(refPhCyc2)) '!!!']); 
                        outON_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2}=pulse_seq(Npts,sw,Bfield,lw,taus,sys,edit_angle_on,...
                            editRFon,editTp,editPhCyc1(EP1),editPhCyc2(EP2),...
                            refRF,refTp,Gx,Gy,x(X),y(Y),refPhCyc1(RP1),refPhCyc2(RP2));
                        outOFF_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2}=pulse_seq(Npts,sw,Bfield,lw,taus,sys,edit_angle_off,...
                            editRFoff,editTp,editPhCyc1(EP1),editPhCyc2(EP2),...
                            refRF,refTp,Gx,Gy,x(X),y(Y),refPhCyc1(RP1),refPhCyc2(RP2));
                    
                        if RP1==1 && RP2==1
                            outON_posxy_epc{X}{Y}{EP1}{EP2}=outON_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2};
                            outOFF_posxy_epc{X}{Y}{EP1}{EP2}=outOFF_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2};
                        else
                            outON_posxy_epc{X}{Y}{EP1}{EP2}=op_addScans(outON_posxy_epc{X}{Y}{EP1}{EP2},outON_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2},xor(RP1==length(refPhCyc1),RP2==length(refPhCyc2)));
                            outOFF_posxy_epc{X}{Y}{EP1}{EP2}=op_addScans(outOFF_posxy_epc{X}{Y}{EP1}{EP2},outOFF_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2},xor(RP1==length(refPhCyc1),RP2==length(refPhCyc2)));
                        end
                    end %end of 1st refocusing phase cycle loop
                end %end of 2nd refocusing phase cycle loop.
                
                if EP1==1 && EP2==1
                    outON_posxy{X}{Y}=outON_posxy_epc{X}{Y}{EP1}{EP2};
                    outOFF_posxy{X}{Y}=outOFF_posxy_epc{X}{Y}{EP1}{EP2};
                else
                    outON_posxy{X}{Y}=op_addScans(outON_posxy{X}{Y},outON_posxy_epc{X}{Y}{EP1}{EP2});
                    outOFF_posxy{X}{Y}=op_addScans(outOFF_posxy{X}{Y},outOFF_posxy_epc{X}{Y}{EP1}{EP2});
                end
                outON_posxy{X}{Y}.x = x;
                outOFF_posxy{X}{Y}.x = x;
                
            end %end of 1st editing phase cycle loop.
        end %end of 2nd editing phase cycle loop.
        
        
        outON=op_addScans(outON,outON_posxy{X}{Y});
        outOFF=op_addScans(outOFF,outOFF_posxy{X}{Y});

    end %end of spatial loop (parfor) in y direction.
end %end of spatial loop (parfor) in x direction.

spect{1} = outOFF_posxy;
spect{2} = outON_posxy;
end