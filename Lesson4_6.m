f = 77e9;
c = 3e8;
lambda = c / f;
phaseIncrement = 45;
antennaSpacing = lambda / 2;

theta = asin(phaseIncrement / 360 * lambda / antennaSpacing);

disp(theta/pi*180);