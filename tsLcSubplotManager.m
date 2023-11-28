classdef tsLcSubplotManager < handle
  properties (Access = public)
    N  % number of rows
    M  % number of cols
    Min = [0.03 0.03];
    Max = [0.97 0.97];
    Gap = [0.03 0.03];
    CellXSize = 300;
    CellYSize = 300;
  end
  properties (Access = private)
    axesMap
    subSubplotManagers
  end
  methods
    
    function obj = tsLcSubplotManager(N, M, varargin)
      obj.N = N;
      obj.M = M;
      args.Min = obj.Min;
      args.Max = obj.Max;
      args.Gap = obj.Gap;
      args.CellXSize = obj.CellXSize;
      args.CellYSize = obj.CellXSize;
      args = easyParseNamedArgs(varargin, args);
      obj.Min = args.Min;
      obj.Max = args.Max;
      obj.Gap = args.Gap;
      obj.CellXSize = args.CellXSize;
      obj.CellYSize = args.CellYSize;
      obj.axesMap = containers.Map('keytype', 'char', 'valuetype', 'any');
      obj.subSubplotManagers = containers.Map('keytype', 'char', 'valuetype', 'any');
    end
    
    function fg = initFigure(obj)
      xSize = obj.CellXSize*obj.M + (obj.Gap(1)*(obj.M - 1));
      ySize = obj.CellYSize*obj.N + (obj.Gap(2)*(obj.N - 1));
      
      fg = figure('position', 20 + [0 0 xSize ySize]);
    end
    
    function ax = createAxes(obj, AxeName, Row, Col, NRow, NCol, varargin)
      if obj.axesMap.isKey(AxeName)
        error(['axes ' AxeName ' already exists']);
      end
      
      args.Gap = obj.Gap;
      args = easyParseNamedArgs(varargin, args);
      Gp = args.Gap;
      
      Xmin   = obj.Min(1);
      Ymin   = obj.Min(2);
      Xgap   = Gp(1);
      Ygap   = Gp(2);
      Xmax   = obj.Max(1) + Xgap;
      Ymax   = obj.Max(2) + Ygap;


      Xsize  = (Xmax - Xmin)./obj.M;
      Ysize  = (Ymax - Ymin)./obj.N;

      Xbox   = Xsize*NCol - Xgap;
      Ybox   = Ysize*NRow - Ygap;

      Xstart = Xmin + Xsize.*(Col - 1);
      Ystart = Ymax - Ysize.*Row;
      
      ax = axes('position',[Xstart,Ystart,Xbox,Ybox]);
      
      obj.axesMap(AxeName) = ax;
    end
    
    function subSubPlotManager = CreateSubSubplotManager(obj, SubSubplotManagerName, Row, Col, NRow, NCol, NSub, MSub, varargin)
      if obj.subSubplotManagers.isKey(SubSubplotManagerName)
        error(['axes ' SubSubplotManagerName ' already exists']);
      end
      
      args.Gap = obj.Gap;
      args.subGap = obj.Gap;
      args = easyParseNamedArgs(varargin, args);
      Gp = args.Gap;
      subGap = args.subGap;
      
      Xmin   = obj.Min(1);
      Ymin   = obj.Min(2);
      Xgap   = Gp(1);
      Ygap   = Gp(2);
      Xmax   = obj.Max(1) + Xgap;
      Ymax   = obj.Max(2) + Ygap;

      Xsize  = (Xmax - Xmin)./obj.M;
      Ysize  = (Ymax - Ymin)./obj.N;

      Xbox   = Xsize*NCol - Xgap;
      Ybox   = Ysize*NRow - Ygap;

      Xstart = Xmin + Xsize.*(Col - 1);
      Ystart = Ymax - Ysize.*Row;
      Xend = Xstart + Xbox;
      Yend = Ystart + Ybox;
            
      subSubPlotManager = tsLcSubplotManager(NSub, MSub, 'Min', [Xstart, Ystart], 'Max', [Xend, Yend], 'Gap', subGap, varargin{:});
      
      obj.subSubplotManagers(SubSubplotManagerName) = subSubPlotManager;
    end
    
    function ax = getAxes(obj, AxeName)
      ax = obj.axesMap(AxeName);
    end
    
    function sm = getSubSubplotManager(obj, AxeName)
      sm = obj.subSubplotManagers(AxeName);
    end
    
    function clear(obj)
      obj.axesMap = containers.Map('keytype', 'char', 'valuetype', 'any');
    end
    
  end
end