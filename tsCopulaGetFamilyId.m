function familyId = tsCopulaGetFamilyId(copulaFamily)

  if strcmpi(copulaFamily, 'gaussian')
    familyId = 1;
  elseif strcmpi(copulaFamily, 't')
    familyId = 2;
  elseif strcmpi(copulaFamily, 'gumbel') 
    familyId = 3;
  elseif strcmpi(copulaFamily, 'clayton') 
    familyId = 4;
  elseif strcmpi(copulaFamily, 'frank')
    familyId = 5;
  else
    error(['copulaFamily not supported: ' copulaFamily]);
  end

end