function handles = tsCopulaYearExtrPlotJcdfTrivar( retLev, jcdf, varargin )
 
args.figPosition = [50, 50, 800, 800];
args.xlbl = 'X';
args.ylbl = 'Y';
args.zlbl = 'Z';
args.fontSize = 12;
args.probRange = [];
args.azimuth = -13.9;
args.elevation = -15.6;
args = tsEasyParseNamedArgs(varargin, args);
figPosition = args.figPosition;
xlbl = args.xlbl;
ylbl = args.ylbl;
zlbl = args.zlbl;
fontSize = args.fontSize;
probRange = args.probRange;
azimuth = args.azimuth;
elevation = args.elevation;

if (size(retLev, 2) ~= 3) || (length(size(jcdf)) ~= 3)
  error('tsCopulaYearExtrPlotSctrTrivar: retLev must be an Nx3 array, jpdf must be XxYxZ');
end

fig = figure('position', figPosition, 'color', 'w');

x = retLev(:,1);
y = retLev(:,2);
z = retLev(:,3);
[yy, xx, zz] = ndgrid(y, x, z);

xslice = [median(x), max(x)];
yslice = max(y);
zslice = [median(z), max(z)];

if ~isempty(probRange)
  jcdf(jcdf < min(probRange)) = min(probRange);
  jcdf(jcdf > max(probRange)) = max(probRange);
end

hslc = slice(xx, yy, zz, jcdf, xslice, yslice, zslice);
xlabel(xlbl);
ylabel(ylbl);
zlabel(zlbl);
view(azimuth, elevation);

ax = gca;
ax.FontSize = fontSize;

cb = colorbar;
cb.FontSize = fontSize;
ylabel(cb, 'Probability', 'fontsize', fontSize);

handles = [fig; hslc(:); cb];

end

