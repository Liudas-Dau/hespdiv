a <- 1
b <- 2
f1 <- function(){
  c <- 3
  d <- 4
  f2 <- function(P){
    assign("calls", calls+1, inherits=TRUE)
    print(calls)
    return(P+c+d)
  }
  calls <- 0
  v <- vector()
  for(i in 1:10){
    v[i] <- f2(P=0)
    c <- c+1
    d <- d+1
  }
  return(v)
}
f1()
?assign

?last
test.fun<-function(data){
  puzzle<-list(data)
  info<-list(length(data))
  call<-0
  sele1<-list()
  sele2<-list()
  environment(f1) <- environment()
  f1<-function(sel){
  call<<-c(call,call[length(call)]+1)
  sel1<-sample(1:length(data),1)
  sel2<-sample(1:length(data),1)
  sele1<<-c(sele1,list(sel1))
  sele2<<-c(sele2,list(sel2))
  if (sel1<=length(sel)){
  info<<-c(info,list(sel1))
  puzzle<<-c(puzzle,list(sel[1:sel1]))
  f1(puzzle[[length(puzzle)]])
  }
  if (sel2<=(length(sel)-sel1)){
    info<<-c(info,list(sel2))
    puzzle<<-c(puzzle,list(sel[(sel1+1):length(sel)][1:sel2]))
    f1(puzzle[[length(puzzle)]])
  }
  }
  f1(data)
  return(list(puzzle,info,call,sele1,sele2))
  }
call<-10

tt<-function() print(environment())
tt()
#basic c veikia geriau

b<-c(a,list(1:5000000000000))

test.fun(1:100)
