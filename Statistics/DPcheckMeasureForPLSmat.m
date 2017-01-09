function x = DPcheckMeasureForPLSmat(x,minX,maxX,filename,Group,Cond,Measure,iG,iC,iS,iF)
ind = find(isnan(x));
if ~isempty(ind)
    warning('NaN values in PLS matrix converted to 0s')
    iG
    Group
    iC
    Cond
    iS
    iF
    filename
    Measure
    x(ind) = 0;
end

ind = find(x<minX);
if ~isempty(ind)
    warning('values below the minimum, converted to minimum')
    iG
    Group
    iC
    Cond
    iS
    iF
    filename
    Measure
    x(ind) = minX;
end

ind = find(x>maxX);
if ~isempty(ind)
    warning('values above the maximum, converted to maximum')
    iG
    Group
    iC
    Cond
    iS
    iF
    filename
    Measure
    x(ind) = maxX;
end