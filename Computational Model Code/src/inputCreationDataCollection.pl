use strict;
use warnings;


my $seed = 1000;
my $evolSeed = 10;
for (my $i=1; $i <= $seed; $i++) {
  for (my $j=1; $j <= $evolSeed; $j++) {
    my $msg = "julia measureDataCollectionForCluster.jl ../dataOutputForPaper/Pop_Seed_".$i."/Pop_Seed_".$i."EvolSeed_".$j;
    print ($msg,"\n");
  }
}
