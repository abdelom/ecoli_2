use strict;
use warnings;


my $infile=$ARGV[0];
my $sufix1=$ARGV[1];
my $sufix2=$ARGV[2];
$infile=~/(.+).sam/;    
my $prefix=$1;
#my $prefix=$ARGV[];
print "$infile\t$prefix\n";
my $outfile1=$prefix."_$sufix1.sam";
my $outfile2=$prefix."_$sufix2.sam";
my $outfile3=$prefix."_$sufix1"."_$sufix2.sam";
my $outfile4=$prefix."_unmapped.sam";
my $outfile5=$prefix."_Plasmid.sam";



my $line;
my $CMI;
my $O1;
my $O2;
my $O3;
my $O4;
my $O5;

my @b;
open($CMI,'<',$infile) or die("open: $!");
open($O1,'>',$outfile1) or die("open: $!");
open($O2,'>',$outfile2) or die("open: $!");
open($O3,'>',$outfile3) or die("open: $!");
open($O4,'>',$outfile4) or die("open: $!");
open($O5,'>',$outfile5) or die("open: $!");
#@SQ    SN:NC_009800.1    LN:4643538
#@SQ    SN:NC_008253.1    LN:4938920
#@PG    ID:bwa    PN:bwa    VN:0.7.17-r1188    CL:bwa mem /media/gluster_disk/iame/U738/mohamed.ghalayini/analyse_polymorphisme_souris_strepto_05112019/HS_536_fasta1.fasta output_cat/60J201_S77_R2_cat.gz output_cat/60J201_S77_R1_cat.gz
#NB501159:29:HH2F2AFXY:2:11101:17496:1035    99    NC_008253.1    4762270    0    151M    =    4762348    229    GATTTGTGCAAGGCGAATAACCGGTGCAACAATCTGACCTGCAAGCATATTAAAAGCAATTAACTGACCAATACTTAAATCCCCGGAAATAACCAGGTGTGCTCCCAACCACAGGTTGATGATCATAAAAGTCTTTTGTATTAAATGTATT    AAAAAE/EEAEEEEEAEAEEEEEEEEAEAEEEEEEA/EEEEEEEEEEEEEEEEEEEEEE/EEEEEAE/EEEAEEEAAE//6AEA<AE/EEEE/EEEEEAEE/EEA66/E6EAEEE6EAEEEEEEAEEA//AE/EEEE//AEE/E/EAEAEA    NM:i:2    MD:Z:128C15C6    MC:Z:151M    AS:i:141    XS:i:141

my $target=0;
my $score;
my $score_alternative;
my $secondmatch;
my $nb1=0;
my $nb2=0;
my $nb3=0;
my $nbp=0;
my $unmapped=0;
my $counter=0;
my  $secondpos;
my $multiplehits;

while(defined ($line = <$CMI>))
	{
        
    if ($counter % 1000000==0)
    {print "$counter\n";}
    $counter++;
	chomp $line;
    if (substr($line,0,1)eq '@')
        {
        next;
        }
    $secondpos=-1;
	@b=split(/\t/,$line);
	#print "$line\n";
        
        #read score
        $line=~/AS:i:(\d+)/;
        $score=$1;
        
        #read score of alternative mapping
        $line=~/XS:i:(\d+)/;
        $score_alternative=$1;
        
        if ($b[2] eq "*" || $score==0)
        {
            $unmapped=$unmapped+1;
            $target=0;
            print $O4 "$line\n";
            next;
        }
        if ($b[2] eq "$sufix1")
            {
            $target=1;
            }
        if ($b[2] eq "$sufix2")
            {
            $target=2;
            }
        if ($b[2] eq "Plasmid")
            {
            $target=4;
            }
        
        #if there is a lesser mapping but relavant enough to be noted we have
        if ($line=~/XA:Z:(.+),.(\d+),.*,\d+;/ && $score_alternative>$score-30)#map genome 
        {
        $secondmatch=$1;
            if ($secondmatch ne $b[2]) #if it matches on the other genome store the ratio
            {
             $secondpos=$2;
            }
            
        #if there are more than two alternative hits
        $multiplehits=0;
            if ($line=~/XA:Z:.+,.\d+,.*,\d+;.+,.\d+,.*,\d+;/)
                {
                    $multiplehits=1;
                    $target=3;# poubelle 
                }
            
        
        }
        
        
        #if the alternative mapping is as good
        if ($score_alternative==$score && $score>0)
        {
            #get the position of the second mapping.
            #if ($line=~/XA:Z:(.+),.\d+,.*,\d+;/)
            #{
            #$secondmatch=$1;
            
                #if ($secondmatch ne $b[2] ) #if it matches on the other genome put it in the common pool
                #{
             #    $target=3;
                #}
            #}
            
            $target=3;
        }
        
        
       
        
        if ($target==1)
        {
            print $O1 "$line\n";#"$b[0]\n$b[9]\n+\n$b[10]\n";
        #print $O1 "$line\n";#"$b[0]\n$b[9]\n+\n$b[10]\n";
            $nb1++; #
        }
        if ($target==2)
        {
            print $O2 "$b[3]\t$secondpos\n";
            #print $O2 "$line\n";#"$b[0]\n$b[9]\n+\n$b[10]\n";
            $nb2++;
        }
        if ($target==3)
        {print $O3 "$line\n";
            $nb3++;
        }
        if ($target==4)
        {print $O5 "$line\n";
            $nbp++;
        }
        
	}
if($secondpos != -1){
print "unmapped\t$sufix1\t$sufix2\tboth\tPlasmid\n";}
print "$unmapped\t$nb1\t$nb2\t$nb3\t$nbp\n";
close $CMI;
close $O1;
close $O2;
close $O3;
close $O4;
close $O5;

exit;

open($O3,'<',$outfile3) or die("open: $!");
open($O1,'>>',$outfile1) or die("open: $!");
open($O2,'>>',$outfile2) or die("open: $!");

#find proba to dsitribute the reads.
my $proba=$nb1/($nb1+$nb2);

$nb1=0;
$nb2=0;
while(defined ($line = <$O3>))
{
chomp $line;
@b=split(/\t/,$line);
if (rand()<$proba)
{
  print $O1 "$b[0]\n$b[9]\n+\n$b[10]\n";
    $nb1++;
}
else
{
print $O2 "$b[0]\n$b[9]\n+\n$b[10]\n";
    $nb2++;
}
}
print "$sufix1\t$sufix2\n";
print "$nb1\t$nb2\n";

system("rm $outfile3\n");
close $O1;
close $O2;
close $O3;

# mutation près