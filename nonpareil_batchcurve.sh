library (Nonpareil)
samples <- read.table('samples.txt', sep='\t', h=T);
attach(samples);
np <- Nonpareil.curve.batch(File, r=R, g=G, b=B, libnames=Name, modelOnly=TRUE);
Nonpareil.legend('bottomright');
detach(samples);
