function copulaFamily = tsCopulaGetFamilyFromId(familyId)

  if familyId == 1
    copulaFamily = 'gaussian';
  elseif familyId == 2
    copulaFamily = 't';
  elseif familyId == 3
    copulaFamily = 'gumbel';
  elseif familyId == 4
    copulaFamily = 'clayton';
  elseif familyId == 5
    copulaFamily = 'frank';
  else
    error(['copula familyId not supported: ' num2str(familyId)]);
  end

end