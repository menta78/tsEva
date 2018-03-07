function cleaned_series = tsRemoveConstantSubseries( srs, stackedValuesCount )

cleaned_series = srs;
[tmp1,tmp2] = tsSameValuesSegmentation(diff(srs),0);
for i=1:length(tmp2);
    ii=tmp2{i};
    if length(ii) >= stackedValuesCount; 
        cleaned_series(ii(2:end))=nan;
    end;
end

end

