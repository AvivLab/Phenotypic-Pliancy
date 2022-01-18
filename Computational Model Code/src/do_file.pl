use strict;
use warnings;


my $seed = 1000;
for (my $i=1; $i <= $seed; $i++) {
    my $msg = "julia codeForDatasetGenerationForPaper.jl $i";
    print ($msg,"\n");
}
