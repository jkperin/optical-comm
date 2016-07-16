function FS=getAgilentWaveform(Agilent)
Acquisition_number = 1; % Define number Waveform Matlab will trigger and import
Memory = num2str(Agilent.WfmLength);
%instrreset();
v=visa('agilent', ['TCPIP0::' Agilent.IPAddr '::inst0::INSTR']);
%% 	Set the scope’s input buffer size 
buffersize=str2double(Memory)*2;
% disp(buffersize);
set (v, 'InputBufferSize', buffersize);
%% Open the instrument with fopen
% disp('Opening Oscilloscope Communication');
fopen(v);
v.Timeout=10;
%% Verify you are talking to the instrument
% Now that v is open,  use the query command with *IDN? to verify you
% are talking to the scope.  You will get a print out of the instrument’s 
% identity
scopeid=query(v,'*IDN?');
% disp(['Oscillosocpe identified: ',scopeid]);
status=get(v,'status');
% disp(status);

fprintf(v,':RUN');
fprintf(v,[':ACQuire:POINts:ANALOG ',Memory]);
fprintf(v,[':ACQuire:INTERPOLATE OFF']);

%Loop acquisitions
for i=0:1:Acquisition_number-1

    scopeOK = str2num(query(v,'ADER?'));
% disp(scopeOK);

%fprintf(v,':DIGITIZE CHANNEL1');

fprintf(v,':DIG');
fprintf(v,':CHANNEL1:DISPLAY ON');
fprintf(v,':CHANNEL2:DISPLAY ON');
fprintf(v,':CHANNEL3:DISPLAY ON');
fprintf(v,':CHANNEL4:DISPLAY ON');

% Need to Wait for Acquisition to complete
scopeOK=0;
while(scopeOK==0);
    scopeOK = str2num(query(v,'ADER?'));
    
end;


%% Get CHANNEL 1 Scaling INfos
% Set variables xinc and yinc to X and Y increment and xor, yor to X and 
% Y origin
% For example:
% yinc=str2double(query(v,':WAVEFORM:YINC?'));
fprintf(v,':WAVEFORM:SOURCE CHAN1');
fprintf(v,':WAVEFORM:FORMAT WORD');
fprintf(v,':WAVEFORM:BYTEORDER LSBFIRST');
fprintf(v,':SYSTEM:HEADER Off');
fprintf(v,':WAVEFORM:STREAMING OFF');
yinc=str2double(query(v,':WAVEFORM:YINC?'));
yor=str2double(query(v,':WAVEFORM:YOR?'));
xinc=str2double(query(v,':WAVEFORM:XINC?'));
xor=str2double(query(v,':WAVEFORM:XOR?'));
fprintf(v,':WAVEFORM:DATA?')
y = binblockread(v,'short');
scopeOK=0;
t=[0:length(y)-1];
%% Scale y with yinc and yor, t with xinc and xorf
y=y*yinc+yor;
t=t*xinc+xor;
%% Create a new figure window (help figure) and plot handle h=plot(t,y)
%figure(1);
%h=plot(t,y);
CH1_Data=[y];
CH1_Time=[t'];
% End of CH1 Data Processing


%% Get CHANNEL 2 Scaling INfos
% Set variables xinc and yinc to X and Y increment and xor, yor to X and 
% Y origin
% For example:
% yinc=str2double(query(v,':WAVEFORM:YINC?'));
fprintf(v,':WAVEFORM:SOURCE CHAN2');
fprintf(v,':WAVEFORM:FORMAT WORD');
fprintf(v,':WAVEFORM:BYTEORDER LSBFIRST');
fprintf(v,':SYSTEM:HEADER Off');
fprintf(v,':WAVEFORM:STREAMING OFF');
yinc=str2double(query(v,':WAVEFORM:YINC?'));
yor=str2double(query(v,':WAVEFORM:YOR?'));
xinc=str2double(query(v,':WAVEFORM:XINC?'));
xor=str2double(query(v,':WAVEFORM:XOR?'));
%% Read the first block of data
% ':WAVEFORM:FORMAT BYTE'
% ':ACQUIRE:TYPE NORM'
% ':WAVEFORM:DATA?'
% Then, use binblockread (use help binblockread if necessary) to read the 
% data into a vector y. Note that the
% precision on binblockread must be int8 to be compatible with
% ':WAVEFORM:FORMAT BYTE'


% disp(['Importing CH2 Waveform N°',num2str(i+1)]);


fprintf(v,':WAVEFORM:DATA?')
y = binblockread(v,'short');
%fread(v);



%b=fread(v);
%endchar=fread(v); % read # Char at End of Binary Block
%disp('Importing command sent');
% Reading Waveform as Binary array of Signed BYTES
% correct type is INT8

scopeOK=0;
%while(scopeOK==0);
 %   scopeOK = str2num(query(v,'*OPC?'));
  %  disp('Checking *OPC?');
   % disp(scopeOK);
    
%end;

%% Create a time vector t which runs from 0 to length(y)-1
t=[0:length(y)-1];
%% Scale y with yinc and yor, t with xinc and xorf
y=y*yinc+yor;
t=t*xinc+xor;
%% Create a new figure window (help figure) and plot handle h=plot(t,y)
%figure(2);
%h=plot(t,y);
CH2_Data=[y];
CH2_Time=[t'];
% End of CH1 Data Processing


%% Get CHANNEL 3 Scaling INfos
% Set variables xinc and yinc to X and Y increment and xor, yor to X and 
% Y origin
% For example:
% yinc=str2double(query(v,':WAVEFORM:YINC?'));
fprintf(v,':WAVEFORM:SOURCE CHAN3');
fprintf(v,':WAVEFORM:FORMAT WORD');
fprintf(v,':WAVEFORM:BYTEORDER LSBFIRST');
fprintf(v,':SYSTEM:HEADER Off');
fprintf(v,':WAVEFORM:STREAMING OFF');
yinc=str2double(query(v,':WAVEFORM:YINC?'));
yor=str2double(query(v,':WAVEFORM:YOR?'));
xinc=str2double(query(v,':WAVEFORM:XINC?'));
xor=str2double(query(v,':WAVEFORM:XOR?'));
%% Read the first block of data
% ':WAVEFORM:FORMAT BYTE'
% ':ACQUIRE:TYPE NORM'
% ':WAVEFORM:DATA?'
% Then, use binblockread (use help binblockread if necessary) to read the 
% data into a vector y. Note that the
% precision on binblockread must be int8 to be compatible with
% ':WAVEFORM:FORMAT BYTE'


% disp(['Importing CH3 Waveform N°',num2str(i+1)]);


fprintf(v,':WAVEFORM:DATA?')
y = binblockread(v,'short');
%fread(v);



%b=fread(v);
%endchar=fread(v); % read # Char at End of Binary Block
%disp('Importing command sent');
% Reading Waveform as Binary array of Signed BYTES
% correct type is INT8

scopeOK=0;
%while(scopeOK==0);
 %   scopeOK = str2num(query(v,'*OPC?'));
  %  disp('Checking *OPC?');
   % disp(scopeOK);
    
%end;

%% Create a time vector t which runs from 0 to length(y)-1
t=[0:length(y)-1];
%% Scale y with yinc and yor, t with xinc and xorf
y=y*yinc+yor;
t=t*xinc+xor;
%% Create a new figure window (help figure) and plot handle h=plot(t,y)
%figure(3);
%h=plot(t,y);
CH3_Data=[y];
CH3_Time=[t'];
% End of CH1 Data Processing


%% Get CHANNEL 4 Scaling INfos
% Set variables xinc and yinc to X and Y increment and xor, yor to X and 
% Y origin
% For example:
% yinc=str2double(query(v,':WAVEFORM:YINC?'));
fprintf(v,':WAVEFORM:SOURCE CHAN4');
fprintf(v,':WAVEFORM:FORMAT WORD');
fprintf(v,':WAVEFORM:BYTEORDER LSBFIRST');
fprintf(v,':SYSTEM:HEADER Off');
fprintf(v,':WAVEFORM:STREAMING OFF');
yinc=str2double(query(v,':WAVEFORM:YINC?'));
yor=str2double(query(v,':WAVEFORM:YOR?'));
xinc=str2double(query(v,':WAVEFORM:XINC?'));
xor=str2double(query(v,':WAVEFORM:XOR?'));
%% Read the first block of data
% ':WAVEFORM:FORMAT BYTE'
% ':ACQUIRE:TYPE NORM'
% ':WAVEFORM:DATA?'
% Then, use binblockread (use help binblockread if necessary) to read the 
% data into a vector y. Note that the
% precision on binblockread must be int8 to be compatible with
% ':WAVEFORM:FORMAT BYTE'


% disp(['Importing CH4 Waveform N°',num2str(i+1)]);


fprintf(v,':WAVEFORM:DATA?')
y = binblockread(v,'short');
%fread(v);



%b=fread(v);
%endchar=fread(v); % read # Char at End of Binary Block
%disp('Importing command sent');
% Reading Waveform as Binary array of Signed BYTES
% correct type is INT8

scopeOK=0;
%while(scopeOK==0);
 %   scopeOK = str2num(query(v,'*OPC?'));
  %  disp('Checking *OPC?');
   % disp(scopeOK);
    
%end;

%% Create a time vector t which runs from 0 to length(y)-1
t=[0:length(y)-1];
%% Scale y with yinc and yor, t with xinc and xorf
y=y*yinc+yor;
t=t*xinc+xor;
%% Create a new figure window (help figure) and plot handle h=plot(t,y)
%figure(4);
%h=plot(t,y);
CH4_Data=[y];
CH4_Time=[t'];
% End of CH1 Data Processing



end;
%END OF LOOP

% disp('Closing Oscilloscope Communication')
fclose(v);
FS = [CH1_Data CH2_Data CH3_Data CH4_Data];
