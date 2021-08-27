my $header;
while (<VCF>) {
        chomp;
        if ($_=~/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t(.+)/) {
                @samples=split(/\s/,$1);
                last;
        }
}
my $nsam=$#samples+1;
my @loci=();
my %lgts={};
my @gts=();
while (<VCF>) {
        chomp;
        my $line=$_;
        if ($line=~/(\S+)\t(\S+)\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(.+)/ ){
                my $locus=$1."_".$2;
                push @loci, $locus;
                @gts=split(/\s/,$3);
                for($g=0; $gt=$gts[$g];$g++) {
                        ($gtt,@rest)=split(/:/,$gt);
                        ($gt1,$gt2)=split(/[|\/]/,$gtt);
#print "$g  $gt  $gt1-$gt2  @samples[$g]\n";
                        if ($gt1 eq ".") {
                                $gt1=0;
                                $gt2=0;
                        }
                        else {
                                $gt1=$gt1+1;
                                $gt2=$gt2+1;
                        }
#print "\t\t$gt1-$gt2  @samples[$g]\n";
                        $lgts{$locus}{$samples[$g]}="0".$gt1."0".$gt2;
                }
        }
}

@samples=sort(@samples);

print "converted from $vcf\n";
print join(", ",@loci), "\n";
my $newpop;
my $oldpop;

foreach $sa (@samples){
        foreach $p (@pops) {
                if ($sa=~/$p/) {$newpop=$p;}
        }
        if ($newpop ne $oldpop) {
                print "POP\n";
                $oldpop=$newpop;
        }
        print $sa, ",";
        foreach $l (@loci) {
                print " ", $lgts{$l}{$sa};
        }
        print "\n";
}
