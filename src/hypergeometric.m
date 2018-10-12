function hyperpval = hypergeometric(N,M,K,x)

hyperpval = nchoosek(M,x)*nchoosek(N-M,K-x)/nchoosek(N,K);
