function dt = tsEvaGetTimeStep( times )
  df = diff(times);
  dt = min(df);
  if dt == 0
      df = df(df ~= 0);
      dt = min(df);
  end
end

