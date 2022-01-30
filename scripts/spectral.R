library(tidyverse)
DATA_PATH = 'data/nbhd_coord_schurch_2020/'
df = read_csv(paste0(DATA_PATH,'CRC_master.csv'))

d = df %>% select(X,Y,ClusterName,spots) %>%
  filter(spots == '15_A') %>%
  select(-spots) %>%
  rename(type = ClusterName)

d$X_scale = (d$X - min(d$X)) / (max(d$X) - min(d$X))
d$Y_scale = (d$Y - min(d$Y)) / (max(d$Y) - min(d$Y))

fft(cbind(d$X_scale,d$Y_scale))

fft0 <- function(z) {
  n <- length(z)
  if(n == 0) return(z)
  k <- 0:(n-1)
  ff <- -2*pi * 1i * k/n
  vapply(1:n, function(h) sum(z * exp(ff*(h-1))), complex(1))
}

dft <- function(z,P,Q) {
  ff <- outer(0:(P-1),-Q:(Q-1),Vectorize(function(p,q) {
    sum(apply(z,1,function(row) {
      exp(-2*pi*1i*(p*row[1]+q*row[2]))
    }))
  }))
  ff
}

z = cbind(X$X_scale,X$Y_scale)

sum(apply(z,1,function(row) {
  exp(-2*pi*1i*(p*row[1]+q*row[2]))
}))

ff <- outer(0:(P-1),-Q:(Q-1),Vectorize(function(p,q) {
  sum(apply(z,1,function(row) {
    exp(-2*pi*1i*(p*row[1]+q*row[2]))
  }))
}))

fft0(cbind(d$X_scale,d$Y_scale))

types = unique(d$type)

X = d %>% filter(type == 'vasculature')
dft(cbind(X$X_scale,X$Y_scale),16,16)
arr = array(0,dim = c(16,32,length(types)))
arr2 = array(0,dim = c(16,32,length(types),length(types)))
for(i in 1:length(types)) {
  t = types[i]
  z = d %>% filter(type == t)
  F_i = dft(cbind(z$X_scale,z$Y_scale),16,16)
  arr[,,i] <- F_i
}
for(i in 1:length(types)) {
  for(j in 1:length(types)) {
    arr2[,,i,j] = arr[,,i]*Conj(arr[,,j])
  }
}
for(i in 1:16) {
  # print(i)
  for(j in 1:32) {
    # print(j)
    # arr2[i,j,,] = solve(arr2[i,j,,])
    dia = diag(arr2[i,j,,])
    arr2[i,j,,] = arr2[i,j,,]^2
    arr2[i,j,,] = sweep(arr2[i,j,,],1,dia,"/")
    arr2[i,j,,] = sweep(arr2[i,j,,],2,dia,"/")
    arr2[i,j,,] = -Mod(sqrt(arr2[i,j,,]))
  }
}

arr2[1,1,,]
