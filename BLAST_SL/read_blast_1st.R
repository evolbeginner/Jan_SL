 e=unlist(read.table("e.list")); p=1-exp(-e); q=p.adjust(p, method="fdr")
