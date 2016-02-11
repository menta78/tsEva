function year = tsYear(timeStamp)

dtvc = datevec(timeStamp);
year = dtvc(:,1);