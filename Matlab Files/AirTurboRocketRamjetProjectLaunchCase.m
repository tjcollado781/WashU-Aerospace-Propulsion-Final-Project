clc
close all
clear all


FanPressure = 1.1:0.1:4
phis = [1 2 10]

SFCs = zeros(length(FanPressure),length(phis));
ISPs = zeros(length(FanPressure),length(phis));
SpecificThrusts = zeros(length(FanPressure),length(phis));
AreaRatios = zeros(length(FanPressure),length(phis));

for i = 1:length(FanPressure)
    for j = 1:length(phis)
        output = AirTurboRocketRamjetIterations(30000,FanPressure(i),phis(j));
        SFCs(i,j) = output(1);
        ISPs(i,j) = output(2);
        SpecificThrusts(i,j) = output(3);
        AreaRatios(i,j) = output(4);
    end
end

plot(FanPressure,SFCs)
title("SFC vs. Fan Pressure Ratio - Launch Case")
xlabel("Fan Pressure Ratio")
ylabel("SFC")
legend('\phi = 1','\phi = 2', '\phi = 10', 'Location', 'NorthEast')

figure
plot(FanPressure,ISPs)
title("Specific Impulse vs. Fan Pressure Ratio - Launch Case")
xlabel("Fan Pressure Ratio")
ylabel("Specific Impulse (s)")
legend('\phi = 1','\phi = 2', '\phi = 10', 'Location', 'NorthWest')

figure
plot(FanPressure,SpecificThrusts)
title("Specific Thrust vs. Fan Pressure Ratio - Launch Case")
xlabel("Fan Pressure Ratio")
ylabel("Specific Thrust")
legend('\phi = 1','\phi = 2', '\phi = 10', 'Location', 'NorthEast')

figure
plot(FanPressure,AreaRatios)
title("Area Ratio vs. Fan Pressure Ratio - Launch Case")
xlabel("Fan Pressure Ratio")
ylabel("Area Ratio")
legend('\phi = 1','\phi = 2', '\phi = 10', 'Location', 'NorthEast')