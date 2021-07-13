function [d,h,k]=readsacfile(sacfile)
% read a single sac file
% a subfunction written by Greenfield & Battenhouse

% Read header first
sacfid=fopen(sacfile,'r');
fparam=fread(sacfid,70,'float32');          %floating values
iparam=fread(sacfid,35,'int32');            %integer values
fseek(sacfid,110*4,-1);
kparam=char(fread(sacfid,24*8,'char')');    %alphanumeric values
[h,k]=assignparam(fparam,iparam,kparam);
for i=1:40
    if h(i)==-12345
        h(i)=NaN;
    end
end;
% Ascertain number of points, preallocate, then read data
np=h(35);
d=zeros(np,1);

% Jump to first data point. Header is 158 words (632 bytes) long.
fseek(sacfid,158*4,-1);  
d(1:np,1)=fread(sacfid,np,'float32'); %make it a column vector
fclose(sacfid);
end

function [h,k]=assignparam(fparam,iparam,kparam)
h = zeros(40,1);
k = cell(4,1);
h(1:3) = fparam(1:3);                   %delta,depmin,depmax,
h(4:5) = fparam(6:7);                   %B,E
h(6) = fparam(8);                       %Event Origin Time
h(7:16) = fparam(11:20);                %T0 - T9
h(17:24) = fparam(32:39);               %STLA,LO,EL,DP,EVLA,LO,EL,DP
h(25:28) = fparam(51:54);               %Dist,Az,BAz,GCArc
h(29:34) = iparam(1:6);                 %Year,Day,Hour,Min,Sec,MillSec
h(35) = iparam(10);                     %Npts
h(36:38) = iparam(16:18);               %IFType,IDep,IZType
h(39:40) = [iparam(23) iparam(25)];     %IevTyp,ISynth
k{1} = deblank(kparam(1:8));            %kstnm
k{2} = deblank(kparam(161:168));        %kcmpnm
k{3} = deblank(kparam(25:32));          %khole
if strcmp(k{3},'-12345')
    k{3} = '';
end
k{4} = deblank(kparam(169:176));        %knetwk
end