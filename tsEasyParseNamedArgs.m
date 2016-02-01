function [argStruct] = tsEasyParseNamedArgs(args, argStruct)
  avlArgs = fieldnames(argStruct);
  for ia = 1:length(avlArgs)
      argName = avlArgs{ia};
      argIndx = find(strcmpi(args, argName));
      if ~isempty(argIndx)
          val = args{argIndx + 1};
          argStruct.(argName) = val;
      end
  end

end

