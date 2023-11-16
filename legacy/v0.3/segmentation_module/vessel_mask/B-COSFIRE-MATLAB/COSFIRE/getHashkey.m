function key = getHashkey(tuple)

primelist = [2 3 5 7 11 13 17 23 29];
key = prod(primelist(1:length(tuple)).^tuple);