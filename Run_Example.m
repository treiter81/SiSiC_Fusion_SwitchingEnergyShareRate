%clear all
%close all

load ExampleData.mat

set(0, 'DefaultFigureRenderer', 'painters');

close all
Energy_ONFUSION(LSopt_Eon_1200A,0) %1st element is the scope data, second a time shift for plotting
Energy_ONFUSION(STD_Eon_1200A,500)
Energy_ONFUSION(STD_Eon_600A,1000)
Energy_ONFUSION(LSopt_Eon_600A,1500)
Energy_ONFUSION(STD_Eon_100A,2000)

xlim([-200 3700])
subplot(211)
legend('V_{ds}','I_{sum}','I_{SLC}')
