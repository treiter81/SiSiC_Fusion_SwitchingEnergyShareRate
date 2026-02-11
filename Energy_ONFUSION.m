function [OUTPUT] = Energy_ONFUSION(input,timeoffset)

%OUTPUT = [Op.Vce Op.Ic Op.E Op.Icpeak Op.Icpeak-Op.Ic Op.dVdt Op.dIdt Op.dIdtSIC Op.dIdtSICandIGBT Op.Enorm share].*[1 1 1e3 1 1 1e-9 1e-9 1e-9 1e-9 1e3 100];

if (nargin<2)
    timeoffset=0;    
end

debug=1;


dat.t = input(:,1);
dat.Vge = input(:,2);
dat.Vce = input(:,3);
dat.Ic = input(:,4);


[l,w] = size(input);
if w<5
    dat.Iigbt = NaN*zeros(l,1);
    dat.Imos = NaN*zeros(l,1);
    dat.Ia = NaN*zeros(l,1);
else
    dat.Iigbt = input(:,5);
    dat.Imos = input(:,6);
    dat.Ia = input(:,7);
end
clear l w


dat.Iigbt = dat.Iigbt + dat.Ia;
dat.Imos = dat.Ic - dat.Iigbt;


param.dvdt_low = 0.1; %typisch 10/90% Werte dvdt
param.dvdt_high = 0.9; %typisch 10/90% Werte dvdt

param.didt_low = 0.1; %typisch 10/90% Werte dvdt
param.didt_high = 0.9; %typisch 10/90% Werte dvdt

param.E_start = 0.1;    %typisch 10% Vge
param.E_stop = 0.02;   %typisch 2% Vce

param.Ipeakfilter = 5; %number of sample points to average

%param.dts = 5e-9;     %time range in s for the Ic derivates (typ 5 to 10ns)
param.fc = 50e6;      %filter the Ic waveform. Recommended to set bandwidth limit of current probe in Hz (typ 30 to 50MHz)
param.dts = 1/param.fc;     %time range 1/bandwidth of current probe

param.peakthreshold = 1e16;   %2nd derivate peak threshold value. Below ignored to cancel noise



%%%%% OFFSET Correction
dat.Vce = dat.Vce - mean(dat.Vce(end-100:end));
dat.Ic = dat.Ic - mean(dat.Ic(1:100)); 
dat.Imos = dat.Imos -  mean(dat.Imos(1:100)); 
dat.Iigbt = dat.Iigbt - mean(dat.Iigbt(1:100)); 

%%%%% Timestep Detection
ts = (dat.t(end) - dat.t(1))/(length(dat.t)-1);
dat.E = ts * cumsum(dat.Vce.*dat.Ic); %[J]
dat.Eigbt = ts * cumsum(dat.Vce.*dat.Iigbt); %[J]
dat.Emos = ts * cumsum(dat.Vce.*dat.Imos); %[J]


%%%%% Operation Point Detection
Op.Vce = mean(dat.Vce(1:100));
[~,posIpeak] = max(dat.Ic);
Op.Icpeak = mean(dat.Ic(posIpeak-int32(param.Ipeakfilter/2):posIpeak+int32(param.Ipeakfilter/2)));

%OP.Ic
[~,posImin] = min(dat.Ic(posIpeak:end)); %find min after current peak
posImin = posImin + posIpeak;
if length(dat.Ic(posImin:end)) > (1e-6/ts)
    Op.Ic = mean(dat.Ic(posImin:(posImin + int32(1e-6/ts))));
else
    Op.Ic = mean(dat.Ic(posImin:end));
end


%dvdt Measurement
pos50V = find(dat.Vce < (Op.Vce/2),1);
poshighV = pos50V - find(dat.Vce(pos50V:-1:1) > (Op.Vce*param.dvdt_high),1);
poslowV  = pos50V + find(dat.Vce(pos50V:1:end) < (Op.Vce*param.dvdt_low),1);
Op.dVdt = (dat.Vce(poshighV) - dat.Vce(poslowV)) / (dat.t(poshighV)-dat.t(poslowV)); %[V/s]


%didt Measurement
pos50I = find(dat.Ic > (Op.Ic/2),1);
poslowI  = pos50I - find(dat.Ic(pos50I:-1:1) < (Op.Ic*param.didt_low),1);
poshighI = pos50I + find(dat.Ic(pos50I:1:end) > (Op.Ic*param.didt_high),1);
poshighIpeak = pos50I + find(dat.Ic(pos50I:1:end) > (Op.Icpeak*param.didt_high),1);
Op.dIdt = (dat.Ic(poshighI) - dat.Ic(poslowI)) / (dat.t(poshighI)-dat.t(poslowI)); %[A/s]


%Check if single or double di/dt slope within continous didt and signal
dat.dIc = [zeros(int32(0.25*param.dts/ts),1); (1/param.dts) .* (dat.Ic(int32(0.5*param.dts/ts)+1:end) - dat.Ic(1:end-int32(0.5*param.dts/ts))); zeros(int32(0.25*param.dts/ts),1)];
% didt = SR * Ic[k+0.5*(SR/fc)] - Ic[k-0.5*(SR/fc)] 
% with SR = 1/Ts ; SR sample rate of Signal
% with k points of signal
% fc bandwidth of current probe 


%filter and make 2nd derivate
fs = 1/ts;
[b,a] = butter(5,param.fc/(fs/2));
dat.dIcfilt = filtfilt(b,a,dat.dIc);
dat.d2Ic = [zeros(int32(0.25*param.dts/ts),1); (1/param.dts) .* (dat.dIcfilt(int32(0.5*param.dts/ts)+1:end) - dat.dIcfilt(1:end-int32(0.5*param.dts/ts))); zeros(int32(0.25*param.dts/ts),1)];

%detect the peaks in 2nd derivate
%[pks,locs,w,p] = findpeaks(dat.d2Ic(1:poshighIpeak),'SortStr','descend');
[pks,locs,w,prom] = findpeaks(dat.d2Ic(1:poshighIpeak));
A.all = [pks,locs,w,prom];
A.prom = sortrows(A.all,4,'descend');           %sort to prominence
A.peaks = sortrows(A.prom(1:2,:),1,'descend');  %sort to peak; only the two highest prominent ones
A.locs = sortrows(A.peaks,2,'ascend');         %sort the two peaks to location

clear pks locs w prom

%pks=A(:,1);
%locs=A(:,2);
%w=A(:,3);
%prom=A(:,4);

%%%%%%%%%%%%check on significance of two slopes
try 
    if (A.peaks(2,1) < param.peakthreshold)
        locs = A.peaks(1,2);
        pks = A.peaks(1,1);
        disp('single slope detected')
    else
        locs = A.locs(:,2);
        pks = A.locs(:,1);        
        disp('double slope detected at current [A]')
        Op.Islc = dat.Ic(locs(2));
        disp(Op.Islc)

    end
catch
    locs = A.peaks(1,2);
    pks = A.peaks(1,1);
    disp('detection was not possible. Continue with single slope')    
end




if (length(locs)==2) && ((1.1*Op.Islc)<(0.9*Op.Icpeak))
    %didt SiC MOSFET Measurement (between initial and the second di/dt
    %slope change 10% Ipeak to 90% Islc
    Op.Islc = dat.Ic(locs(2));
    pos_slc = locs(2);

    %pos50I = find(dat.Ic > (Op.Ic/2),1);    
    poslowI_SIC  = pos_slc - find(dat.Ic(pos_slc:-1:1) < (0.1 * Op.Icpeak),1);
    poshighI_SIC = pos_slc + find(dat.Ic(pos_slc:1:end) > (0.9 * Op.Islc),1);
    Op.dIdtSIC = (dat.Ic(poshighI_SIC) - dat.Ic(poslowI_SIC)) / (dat.t(poshighI_SIC)-dat.t(poslowI_SIC)); %[A/s]
    %dist1 = uint32(0.1*locs(2)-locs(1));
    %Op.dIdtSIC = (dat.Ic(locs(2)-dist1) - dat.Ic(locs(1)+dist1)) / (dat.t(locs(2)-dist1)-dat.t(locs(1)+dist1)); %[A/s]
    
    %didt SiC+IGBT MOSFET Measurement (between the second di/dt and max
    %110 Iscl to 90% Ipeak
    
    poslowI_BOTH  = pos_slc + find(dat.Ic(pos_slc:1:end) > (1.1 * Op.Islc),1);
    poshighI_BOTH = pos_slc + find(dat.Ic(pos_slc:1:end) > (0.9 * Op.Icpeak),1);
    Op.dIdtSICandIGBT = (dat.Ic(poshighI_BOTH) - dat.Ic(poslowI_BOTH)) / (dat.t(poshighI_BOTH)-dat.t(poslowI_BOTH)); %[A/s]
    %dist2 = uint32(0.1*poshighIpeak-locs(2));
    %Op.dIdtSICandIGBT = (dat.Ic(poshighIpeak-dist2) - dat.Ic(locs(2)+dist2)) / (dat.t(poshighIpeak-dist2)-dat.t(locs(2)+dist2)); %[A/s]
    
    %check if di/dt plausible; if not just single value
    if (Op.dIdtSICandIGBT<Op.dIdt)
        Op.dIdtSIC = Op.dIdt;
        Op.dIdtSICandIGBT = NaN;
        locs=locs(1); %clear the other values
        disp('double slope was not plausible and discarded')
    end

else % in case no dual slope only SIC di/dt with IEC nomenclature
    Op.dIdtSIC = Op.dIdt;
    Op.dIdtSICandIGBT = NaN;
end



figure(1)
 
%Enorm Measurement
Op.E = max(dat.E); 
Op.Vge = dat.Vge(end);
Op.Estart = param.E_start*Op.Vge;
Op.Estop = param.E_stop*Op.Vce;
Op.Ei1 = find(dat.Vge > Op.Estart,1);
Op.Ei2  = pos50V + find(dat.Vce(pos50V:1:end) < (Op.Estop),1);
Op.Enorm = dat.E(Op.Ei2) - dat.E(Op.Ei1);


%E Fusion Energy share rates
if (length(locs)==2)
    Op.ESIC1 = dat.E(max(locs)) - dat.E(Op.Ei1); %until slope changing point only SIC is turning on 

    SICI = (Op.dIdtSIC/Op.dIdtSICandIGBT)*mean([dat.Ic(max(locs)) Op.Icpeak]);     %SIC switching the current from slope change until peak with a share rate dependent on di/dt share rates
    IBGTI = (1-(Op.dIdtSIC/Op.dIdtSICandIGBT))*mean([0 Op.Icpeak-dat.Ic(locs(2))]); %IGBT switching the current from 0A to difference Ipeak-slope change with a share rate dependent on di/dt share rates

    Op.ESIC2 = (SICI/(SICI+IBGTI)) * (dat.E(Op.Ei2) - dat.E(max(locs)));           %The energy is now set depending on current share rates
    Op.ESIC = Op.ESIC1+Op.ESIC2;                                                %SIC Energy is from 1st phase only SIC switching until both switching with the calculated share rates

    Op.EIGBT = (1-(SICI/(SICI+IBGTI))) * (dat.E(Op.Ei2) - dat.E(max(locs)));       %IGBT switches only after the slope changing
    
    share = Op.ESIC/Op.Enorm;                                                   %calculates the share rates of SIC turn-on energy related to complete Eon

else
    Op.ESIC = Op.Enorm;    %all to SIC
    share = 1;
end



OUTPUT = [Op.Vce Op.Ic Op.E Op.Icpeak Op.Icpeak-Op.Ic Op.dVdt Op.dIdt Op.dIdtSIC Op.dIdtSICandIGBT Op.Enorm share].*[1 1 1e3 1 1 1e-9 1e-9 1e-9 1e-9 1e3 100];
%disp('EnergySIC')
%disp(Op.ESIC1*1e3)


%Time Offset
dat.t = 1e9*dat.t;  %time in ns
dat.t = dat.t - dat.t(poslowI) + timeoffset;
%dat.t = dat.t - dat.t(find(dat.Vge>1+mean(dat.Vge(1:100)),1));

%Offset for display
dat.E = 1e3*dat.E;  %energy in mJ
dat.E = dat.E - dat.E(Op.Ei1);



ax(1) = subplot(411);
hold on
plot(dat.t,dat.Vge)
ylabel('Vgs [V]')
grid on

ax(2) = subplot(412);
hold on
plot(dat.t,dat.Vce)
ylabel('V_{ds} [V]')
grid on


ax(3) = subplot(413);
hold on
plot(dat.t,dat.Ic)
plot(dat.t(locs),dat.Ic(locs),'r*')
ylabel('I_{sum} [A]')
grid on


ax(4) = subplot(414);
hold on
plot(dat.t,dat.E)
ylabel('E [mJ]')
xlabel('time [ns]')
grid on

linkaxes(ax,'x');
set(gcf,'Position',[10 10 1000 800]);



%debug = 1;
if (debug==1)

    figure(10)
    set(gcf,'Position',[100 100 600 500])
    axd(1) = subplot(211);
    hold on
    plot(dat.t,dat.Ic,'b-','LineWidth',0.75)
    ylabel('I_{sum} [A]')        
    grid on
    ylim([-100 1600])
    limits = axd(1).YLim;

    yyaxis right
    plot(dat.t(1:poshighIpeak),dat.dIc(1:poshighIpeak),'g-')
    plot(dat.t(1:poshighIpeak),dat.dIcfilt(1:poshighIpeak),'r-','LineWidth',1)
    %plot(dat.t,dat.dIc,'g-')
    %plot(dat.t,dat.dIcfilt,'r-','LineWidth',1)
    legend('I_{sum}','di/dt','filtered di/dt','Location','southeast')
    ylabel('di/dt [A/s]')
    xlabel('time [ns]')
    grid on
    axd(1).YLim = limits.*1e7;


    axd(2) = subplot(212);
    hold on
    plot(dat.t,dat.Ic,'b-','LineWidth',0.75)
    axd(2).YLim = limits;
    ylabel('I_{sum} [A]')
    
    yyaxis right
    plot(dat.t(1:poshighIpeak),dat.d2Ic(1:poshighIpeak),'r-')
    plot(dat.t(locs),pks,'r*','MarkerSize',8)

    legend('I_{sum}','d^2i/dt^2','maxima d^2i/dt^2','Location','southeast')
    ylabel('d^2i/dt^2 [A/s^2]')
    xlabel('time [ns]')
    grid on
    axd(2).YLim = limits.*1e14;

    % Create textbox
    annotation(gcf,'textbox',[0.01 0.956 0.044 0.034],'String',{'a)'},'LineStyle','none','FontWeight','bold','FontSize',12,'FitBoxToText','off');
    annotation(gcf,'textbox',[0.01 0.482 0.044 0.034],'String',{'b)'},'LineStyle','none','FontWeight','bold','FontSize',12,'FitBoxToText','off');



    figure(11)
    set(gcf,'Position',[100 100 600 500])
    axd(3) = subplot(211);
    hold on
    plot(dat.t,dat.Vce,'g-')
    plot(dat.t,dat.Ic,'b-')
    

    
    if (length(locs)==2)
        plot(dat.t(locs(2)),dat.Ic(locs(2)),'r*','MarkerSize',8)  
        plot([dat.t(poslowI_SIC), dat.t(poshighI_SIC)],[dat.Ic(poslowI_SIC), dat.Ic(poshighI_SIC)],'m--')  
        plot([dat.t(poslowI_BOTH), dat.t(poshighI_BOTH)],[dat.Ic(poslowI_BOTH), dat.Ic(poshighI_BOTH)],'c--')
        legend('V_{ds}','I_{sum}','I_{SLC}',['di/dt_{SiC} = ',num2str(round(1e-9*Op.dIdtSIC,2)),'kA/us'],['di/dt_{SiC&IGBT} = ',num2str(round(1e-9*Op.dIdtSICandIGBT,2)),'kA/us'],'Location','southeast')        
    else
        plot([dat.t(poslowI), dat.t(poshighI)],[dat.Ic(poslowI), dat.Ic(poshighI)],'m--') 
        legend('V_{ds}','I_{sum}',['di/dt_{SiC} = ',num2str(round(1e-9*Op.dIdtSIC,2)),'kA/us'],'Location','southeast')
    end    
    ylabel('I_{sum} [A], V_{ds} [V]')    
    xlabel('time [ns]')
    grid on

  

    axd(4) = subplot(212);
    hold on
    plot(dat.t(Op.Ei1:Op.Ei2),dat.E(Op.Ei1:Op.Ei2),'b-')
    plot(dat.t(Op.Ei2),dat.E(Op.Ei2),'b*')

    %plot(dat.t(locs),dat.E(locs),'r*','MarkerSize',8)
    plot(dat.t(Op.Ei2),1e3.*Op.ESIC,'k*','LineWidth',1)       
    plot(dat.t(Op.Ei1:Op.Ei2),1e3.*dat.Emos(Op.Ei1:Op.Ei2),'m-')
    plot(dat.t(Op.Ei2),1e3.*dat.Emos(Op.Ei2),'m*')


    if (length(locs)==2)
        plot(dat.t(Op.Ei2),1e3.*Op.EIGBT,'k+','LineWidth',1)
        plot(dat.t(Op.Ei1:Op.Ei2),1e3.*dat.Eigbt(Op.Ei1:Op.Ei2),'c-')    
        plot(dat.t(Op.Ei2),1e3.*dat.Eigbt(Op.Ei2),'c+')    
    end

        

    legend('E_{on,(SiC&IGBT)}','E_{on,(SiC&IGBT)}','Extracted E_{on,(SiC)}','Measured E_{on,(SIC)}','Measured E_{on,(SIC)}','Extracted E_{on,(IGBT)}','Measured E_{on,(IGBT)}','Measured E_{on,(IGBT)}','Location','southeast')
    ylabel('Eon [mJ]')    
    xlabel('time [ns]')
    grid on

    linkaxes(axd,'x')
    axd(1).XLim=[-200 600];

    % Create textbox
    annotation(gcf,'textbox',[0.01 0.956 0.044 0.034],'String',{'a)'},'LineStyle','none','FontWeight','bold','FontSize',12,'FitBoxToText','off');
    annotation(gcf,'textbox',[0.01 0.482 0.044 0.034],'String',{'b)'},'LineStyle','none','FontWeight','bold','FontSize',12,'FitBoxToText','off');
    
end


end

