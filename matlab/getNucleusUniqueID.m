function uniqueID = getNucleusUniqueID(nr,i,k)
    idx = (nr.eggChamberID == i ) & (nr.nucID == k );

    uniqueID = nr.nucUniqueID(idx);

end