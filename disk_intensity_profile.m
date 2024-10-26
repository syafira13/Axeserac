function Ie = disk_intensity_profile(r,param)

yg = param(1); mu = param(2); sig = param(3);

Ie = exp(-((yg + asinh((r-mu)/sig)).^2)/2)./(sqrt(((r-mu).^2)+sig^2));