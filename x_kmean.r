a<-scan("kmeanfea.txt")
x <- matrix(a, ncol = 2000, byrow=TRUE)
set.seed(1357)
km.res=kmeans(x,50,iter.max = 1000, nstart = 25,
       algorithm = c("Hartigan-Wong", "Lloyd", "Forgy",
                     "MacQueen"), trace=FALSE)
print(km.res)
