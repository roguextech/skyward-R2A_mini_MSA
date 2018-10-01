% Auto_Matrices


% This program is used to determine the stability of a rocket through the
% change of some of the parameters of the fins.
% This require in the path for matlab, datcom.exe and datcom_parser.
% The program compiles a datcom input file, generate a datcom output file
% and converts it to a matrix.

%--------------------------------------------------------------------------

% Definition of the INPUT for the code:
% Base : Flag for base plume interaction
% Chord : Panel chord at each semi-span location
% Dnose : Nose diameter at base
% Dcenter : Centerbody diameter at base
% Dexit : Nozzle diameter for base drag calculation
% Lref : Longitudinal reference length
% Latref : Lateral reference length
% Ler : Leading edge radius at each span station
% Lmaxu : Fraction of chord from leading edge to max thickness of upper surface
% Lflatu : Fraction of chord of constant thickness  section of upper surface
% Lnose : Nose length
% Lcenter : Centerbody length
% Npanel : Numbers of panels
% Phif : Angle from each panel
% Sta : Chord station used in measuring sweep (0 = Leading edge)
% Sectype: Cross-section type
% Sspan : Semi-span location
% Sref : Reference area
% Tnose : Nose shape
% XcgF : Longitudinal position of Gravity Center in wet condition
% XcgE : Longitudinal position of Gravity Center in dry condition
% Xle : Distance from missile nose to chord leading edge at each semi-span location
% Zupper : Thickness to chord ratio of upper surface

%--------------------------------------------------------------------------

tic
clear 
close all
clc

%% Flight Conditions

Mach=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1];
Alpha=[-15 -10 -7.5 -5 -2.5 -1.5 -1 -0.5 0.0 0.5 1 1.5 2.5 5 7.5 10 15];
Beta=[-5 0 5];
Alt=[0 500 1000 1500 2000 2500];
Nmach=length(Mach);
Nalpha=length(Alpha);
Nbeta = length(Beta);
Nalt = length(Alt);

%% Design Parameters
% Starting the loop for various dimension of the fins (change these if you
% need some other parameters)
radius= 0.045; % radius of the body
height= [0.075]; 
LSspan=radius+height;
Chord_root=[0.16]; %put the values for different chords
Chord_tip = [0.14];
d = 0; %distance from end of rocket

%% Fixed Parameters
rocket_length = 2.06;

% % Reference Quantities
Xcgf=1.31; 
Xcge=1.15;
% Sref=0.02378;
% Lref=radius*2;
% Latref=radius*2;

% Axisymmetric Body Geometry
Lnose=0.30;
Dnose=0.09;
Lcenter=1.15;
Dcenter=0.09;
Tnose = 'KARMAN'; % nose shape
Dexit = 0.08;

% Define Fin Set n
Npanel=3;
Phif=[0 120 240];
Ler=2*0.0015;
Sta=0; % It means that sweep is measured at the leading edge
Sweep=35;
Sectype = 'ARC'; % Section type (ARC or HEX)

%% Creation of for005.dat

for L1=1:length(LSspan)
    for L2=1:length(Chord_root)
        for L3=1:length(Chord_tip)
            % Create folder for your case
            first_t=height(L1);
            second_t=Chord_root(L2);     
            third_t = Chord_tip(L3);
            Xle=rocket_length - d - Chord_root(L2);
            Sspan=[radius LSspan(L1)];
            h_fin=Sspan(2)-Sspan(1);
            newdir = sprintf('H_%g-C1_%g-C1_%g-D_%g',first_t,second_t,third_t, d);
            mkdir(fullfile(newdir)); 
        
            for Xcg=[Xcgf Xcge]
                % Define Fin Set n (Variables)
                Chord=[Chord_root(L2) Chord_tip(L3)];
                Zupper=[2e-3/Chord(1) 2e-3/Chord(2)];
                Lmaxu=[3e-2/Chord(1) 3e-2/Chord(2)];
                Lflatu=[(Chord(1)-0.06)/Chord(1) (Chord(2)-0.06)/Chord(2)];
                % Creating for005.dat file from previous data
                fid = fopen(strcat(pwd,'\for005.dat'),'w+'); %w+ open or create file for reading and writing; discard existing contents
                % Flight Conditions
                fprintf(fid,'\n $FLTCON\r\n');
                fprintf(fid, '  BETA=0.,\r\n');
                fprintf(fid, '  ALT=%g*0.,\r\n',Nalt);
                fprintf(fid, '  NMACH=%g.',Nmach);
                fprintf(fid, ',');
                fprintf(fid, '\r\n');
                fprintf(fid, '  MACH=');
                aux = min(10,Nmach);
                for i=1:aux
                    fprintf(fid, '%.1f',Mach(i));
                    fprintf(fid, ',');
                end
                if Nmach > 10
                   fprintf(fid, '\r\n');
                   fprintf(fid,  '  MACH(11)=');
                   for i=11:Nmach
                        fprintf(fid, '%.1f',Mach(i));
                        fprintf(fid, ',');
                    end
                end
                
                fprintf(fid, '\r\n');
                fprintf(fid, '  NALPHA=%g.',Nalpha);
                fprintf(fid, ',');
                fprintf(fid, '\r\n');
                fprintf(fid, '  ALPHA=');
                aux = min(10,Nalpha);
                for i=1:aux
                    fprintf(fid, '%.1f',Alpha(i));
                        fprintf(fid, ',');
                 end
                 if Nalpha > 10
                    fprintf(fid, '\r\n');
                    fprintf(fid,  '  ALPHA(11)=');
                    for i=11:Nalpha
                        fprintf(fid, '%.1f',Alpha(i));
                        fprintf(fid, ',');
                    end
                 end
                fprintf(fid, '$');
                
                % Reference Quantities
                fprintf(fid, '\r\n $REFQ\r\n');
                fprintf(fid, '  XCG=');
                fprintf(fid, '%.5f', Xcg);
                % Uncomment this if SREF, LREF, LATREF are specified
%                 fprintf(fid,',\r\n');
%                 fprintf(fid, '  SREF=');
%                 fprintf(fid, '%.5f,\r\n', Sref);
%                 fprintf(fid, '  LREF=');
%                 fprintf(fid, '%.3f,\r\n', Lref);
%                 fprintf(fid, '  LATREF=');
%                 fprintf(fid, '%.3f', Latref);
                fprintf(fid, ',');
                fprintf(fid, '$');

                % Axisymmetric Body Geometry
                fprintf(fid, '\r\n $AXIBOD\r\n');
                fprintf(fid, '  TNOSE=%s',Tnose);
                fprintf(fid, ',');
                fprintf(fid, '\r\n');
                fprintf(fid, '  LNOSE=');
                fprintf(fid, '%.3f,\r\n', Lnose);
                fprintf(fid, '  DNOSE=');
                fprintf(fid, '%.3f,\r\n', Dnose);
                fprintf(fid, '  LCENTR=');
                fprintf(fid, '%.3f,\r\n', Lcenter);
                fprintf(fid, '  DCENTR=');
                fprintf(fid, '%.3f,\r\n', Dcenter);
                fprintf(fid, '  DEXIT=%.3f',Dexit);
                fprintf(fid, ',');
                fprintf(fid, '\r\n');
                fprintf(fid, '  BASE=.FALSE.,$');
                % Finset
                fprintf(fid, '\r\n $FINSET1\r\n');
                fprintf(fid, '  XLE=');
                fprintf(fid, '%.3f,\r\n', Xle);
                fprintf(fid, '  NPANEL=');
                fprintf(fid, '%.1f,\r\n', Npanel);
                fprintf(fid, '  PHIF=');
                for i=1:length(Phif)
                    fprintf(fid, '%.1f', Phif(i));
                    fprintf(fid, ',');
                end
                fprintf(fid, '\r\n');
                fprintf(fid, '  LER=2*');
                fprintf(fid, '%.4f,\r\n', Ler);
                fprintf(fid, '  SWEEP=');
                fprintf(fid, '%.1f,\r\n', Sweep);
                fprintf(fid, '  STA=');
                fprintf(fid, '%.1f,\r\n', Sta);
                fprintf(fid, '  SSPAN=');
                for i=1:length(Sspan)
                    fprintf(fid, '%.3f', Sspan(i));
                    fprintf(fid, ',');
                end
                fprintf(fid, '\r\n');
                fprintf(fid, '  CHORD=');
                for i=1:length(Chord)
                    fprintf(fid, '%.3f', Chord(i));
                    fprintf(fid, ',');
                end
                fprintf(fid, '\r\n');
                fprintf(fid, '  SECTYP=%s,',Sectype);
                fprintf(fid, '\r\n');
                fprintf(fid, '  ZUPPER=');
                for i=1:length(Zupper)
                    fprintf(fid, '%.4f', Zupper(i));
                    fprintf(fid, ',');
                end
                % Uncommen this section if you want HEX type 
%                 fprintf(fid, '\r\n');
%                 fprintf(fid, '  LMAXU=');
%                 for i=1:length(Lmaxu)
%                     fprintf(fid, '%.4f', Lmaxu(i));
%                     fprintf(fid, ',');
%                 end
%                 fprintf(fid, '\r\n');
%                 fprintf(fid, '  LFLATU=');
%                 for i=1:length(Lflatu)
%                     fprintf(fid, '%.4f', Lflatu(i));
%                     fprintf(fid, ',');
%                 end
                fprintf(fid, '$\r\n');

                % Options
                fprintf(fid, 'DERIV RAD\r\n');
                fprintf(fid, 'DIM M\r\n');
                fprintf(fid, 'DAMP\r\n');
                fprintf(fid, 'SAVE\r\n');
                fprintf(fid, 'NEXT CASE\r\n');
                % Cases
                for j=1:Nalt
                    for k=1:Nbeta
                        if Beta(k)==-5 && Alt(j)==0
                        else
                            fprintf(fid,' $FLTCON\r\n');
                            fprintf(fid,'  BETA=');
                            fprintf(fid, '%.1f,\r\n', Beta(k));
                            fprintf(fid,'  ALT=');
                            fprintf(fid, '%.1f', Alt(j));
                            fprintf(fid, '$\r\n');
                            fprintf(fid, 'DERIV RAD\r\n');
                            fprintf(fid, 'DIM M\r\n');
                            fprintf(fid, 'DAMP\r\n');
                            fprintf(fid, 'SAVE\r\n');
                            fprintf(fid, 'NEXT CASE\r\n');
                        end
                    end
                end
                fclose(fid);
                
                %% Creating .dat files+parsing
                disp('Executing datcom.exe');
                dos('datcoming');
                pause(20);
                disp('Creating matrix');
                dos('parsing');
                value=0;
                while value==0
                    value=exist('for006.mat','file');
                    pause(1);
                end
            
                s=sprintf('%s', pwd,'\', newdir);
                if Xcg==Xcgf
                    movefile(strcat(pwd,'\for006.mat'),strcat(s,'\R2A_full.mat'));
                    copyfile(strcat(pwd,'\for005.dat'),strcat(s,'\for005_f.dat'));
                else
                    movefile(strcat(pwd,'\for006.mat'),strcat(s,'\R2A_empty.mat'));
                    copyfile(strcat(pwd,'\for005.dat'),strcat(s,'\for005_e.dat'));
                end
                
            end
        end  
    end
end

beep %gives an alarm when the program stopped working
toc