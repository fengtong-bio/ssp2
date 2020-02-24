#!/usr/bin/perl -w 
#
#Author: Feng Tong
#
my ($WORK_NAME,$WORK_PATH,$KMER_COUNT_PATH,$FILTERX_PATH,$SOAPDENOVO_PATH,$BWA_PATH,$SAMTOOLS_PATH,$THREAD,$G1_SAMPLE_NUMBER,$G2_SAMPLE_NUMBER,$GAP_LENGTH,$MIX_LENGTH,$CONTIG_SPECIFIC_RATIO,$SPECIFIC_KMER_RATIO);
my ($now,$now_path,@G1_1,@G1_2,@G2_1,@G2_2);
my $G1_n = 0;
my $G2_n = 0;
open CFG,"$ARGV[0]" || die;
while (<CFG>){
	chomp;
	$now = $_;
	if ($now =~ /INSTALL_PATH=/){
		$now =~ /INSTALL_PATH=(\S+)/;
		$INSTALL_PATH=$1;
	}

	if ($now =~ /WORK_NAME=/){
		$now =~ /WORK_NAME=(\S+)/;
		$WORK_NAME=$1;
	}
	
	if ($now =~ /WORK_PATH=/){
		$now =~ /WORK_PATH=(\S+)/;
		$WORK_PATH=$1;
	}
	
	if ($now =~ /CFG_PATH=/){
		$now =~ /CFG_PATH=(\S+)/;
		$CFG_PATH=$1;
	}

	if ($now =~ /KMER_COUNT_PATH=/){
		$now =~ /KMER_COUNT_PATH=(\S+)/;
		$KMER_COUNT_PATH=$1;
	}
	
	if ($now =~ /FILTERX_PATH=/){
		$now =~ /FILTERX_PATH=(\S+)/;
		$FILTERX_PATH=$1;
	}
	
	if ($now =~ /SOAPDENOVO_PATH=/){
		$now =~ /SOAPDENOVO_PATH=(\S+)/;
		$SOAPDENOVO_PATH=$1;
	}
	
	if ($now =~ /BWA_PATH=/){
		$now =~ /BWA_PATH=(\S+)/;
		$BWA_PATH=$1;
	}
	
	if ($now =~ /SAMTOOLS_PATH=/){
		$now =~ /SAMTOOLS_PATH=(\S+)/;
		$SAMTOOLS_PATH=$1;
	}
	
	if ($now =~ /THREAD=/){
		$now =~ /THREAD=(\S+)/;
		$THREAD=$1;
	}
	
	if ($now =~ /G1_SAMPLE_NUMBER=/){
		$now =~ /G1_SAMPLE_NUMBER=(\S+)/;
		$G1_SAMPLE_NUMBER=$1;
	}
	
	if ($now =~ /G2_SAMPLE_NUMBER=/){
		$now =~ /G2_SAMPLE_NUMBER=(\S+)/;
		$G2_SAMPLE_NUMBER=$1;
	}
	
	if ($now =~ /GAP_LENGTH=/){
		$now =~ /GAP_LENGTH=(\S+)/;
		$GAP_LENGTH=$1;
	}
	
	if ($now =~ /MIX_LENGTH=/){
		$now =~ /MIX_LENGTH=(\S+)/;
		$MIX_LENGTH=$1;
	}
	
	if ($now =~ /CONTIG_SPECIFIC_RATIO=/){
		$now =~ /CONTIG_SPECIFIC_RATIO=(\S+)/;
		$CONTIG_SPECIFIC_RATIO=$1;
	}
	
	if ($now =~ /SPECIFIC_KMER_RATIO=/){
		$now =~ /SPECIFIC_KMER_RATIO=(\S+)/;
		$SPECIFIC_KMER_RATIO=$1;
	}
	
	if ($now =~ /G1_1_PATH/){
		$now =~ /G1_1_PATH=(\S+)/;
		$now_path=$1;
		$G1_1[$G1_n] = $now_path;
	}
	
	if ($now =~ /G1_2_PATH=/){
		$now =~ /G1_2_PATH=(\S+)/;
		$now_path=$1;
		$G1_2[$G1_n] = $now_path;
		$G1_n++;
	}
	
	if ($now =~ /G2_1_PATH=/){
		$now =~ /G2_1_PATH=(\S+)/;
		$now_path=$1;
		$G2_1[$G2_n] = $now_path;
	}
	
	if ($now =~ /G2_2_PATH=/){
		$now =~ /G2_2_PATH=(\S+)/;
		$now_path =$1;
		$G2_2[$G2_n] = $now_path;
		$G2_n++;
	}

}

my $sit = 0;

my $dir_path = "$WORK_PATH" . "/kmer_count";
if (-e $dir_path){
print STDERR "\nkmer_count exist!\n";
}else{
system("mkdir $dir_path");
}

$dir_path = "$WORK_PATH" . "/kmer_reads";
if (-e $dir_path){
print STDERR "\nkmer_reads exist!\n";
}else{
system("mkdir $dir_path");
}

$dir_path = "$WORK_PATH" . "/G1_SOAPdenovo";
if (-e $dir_path){
print STDERR "\nG1_SOAPdenovo exist!\n";
}else{
system("mkdir $dir_path");
}

$dir_path = "$WORK_PATH" . "/G2_SOAPdenovo";
if (-e $dir_path){
print STDERR "\nG2_SOAPdenovo exist!\n";
}else{
system("mkdir $dir_path");
}

$dir_path = "$WORK_PATH" . "/G1_BWA";
if (-e $dir_path){
print STDERR "\nG1_BWA exist!\n";
}else{
system("mkdir $dir_path");
}

$dir_path = "$WORK_PATH" . "/G2_BWA";
if (-e $dir_path){
print STDERR "\nG2_BWA exist!\n\n";
}else{
system("mkdir $dir_path");
}


print "cd $WORK_PATH\n";
my $filterx = "$FILTERX_PATH -k s -1 cnt>=1";
for (;$sit < $G1_n;$sit++){
print "$KMER_COUNT_PATH -l 21 -f FQ -o $WORK_PATH" . "/kmer_count/G1.$sit.kmer.out -i $G1_1[$sit] -i $G1_2[$sit] -G 10G\n";
$filterx .= " $WORK_PATH" . "/kmer_count/G1.$sit.kmer.out:grp=1";
}
$sit = 0;
for (;$sit < $G2_n;$sit++){
print "$KMER_COUNT_PATH -l 21 -f FQ -o $WORK_PATH" . "/kmer_count/G2.$sit.kmer.out -i $G2_1[$sit] -i $G2_2[$sit] -G 10G\n";
$filterx .= " $WORK_PATH" . "/kmer_count/G2.$sit.kmer.out:grp=1";
}
print "\n\n$filterx > $WORK_PATH" . "/all_kmer.out\n";

print "cd $WORK_PATH/kmer_reads\n";
print "perl $INSTALL_PATH" . "/step1.pl $CFG_PATH $WORK_PATH" . "/all_kmer.out\n";

open G1CFG,">$WORK_PATH" . "/G1.cfg" || die;
open G2CFG,">$WORK_PATH" . "/G2.cfg" || die;
print G1CFG "max_rd_len=150\n[LIB]\navg_ins=300\nreverse_seq=0\nrank=1\nasm_flags=3\nrd_len_cutof=150\npair_num_cutoff=3\nmap_len=32\n";
print G2CFG "max_rd_len=150\n[LIB]\navg_ins=300\nreverse_seq=0\nrank=1\nasm_flags=3\nrd_len_cutof=150\npair_num_cutoff=3\nmap_len=32\n";
$sit = 0;
for (;$sit < $G1_n;$sit++){
print G1CFG "q1=$WORK_PATH" . "/kmer_reads/G1_$sit.R1.fq.gz\nq2=$WORK_PATH" . "/kmer_reads/G1_$sit.R2.fq.gz\n";
}

$sit = 0;
for (;$sit < $G2_n;$sit++){
print G2CFG "q1=$WORK_PATH" . "/kmer_reads/G2_$sit.R1.fq.gz\nq2=$WORK_PATH" . "/kmer_reads/G2_$sit.R2.fq.gz\n";
}

print "cd $WORK_PATH/G1_SOAPdenovo\n";
print "$SOAPDENOVO_PATH all -s ../G1.cfg -o G1_SOAPdenovo -p $THREAD\n";


print "cd $WORK_PATH/G2_SOAPdenovo\n";
print "$SOAPDENOVO_PATH all -s ../G2.cfg -o G2_SOAPdenovo -p $THREAD\n";

print "cd $WORK_PATH\n";
print "ln -s ./G1_SOAPdenovo/G1_SOAPdenovo.scafSeq G1.fa\nln -s ./G2_SOAPdenovo/G2_SOAPdenovo.scafSeq G2.fa\n";

print "$BWA_PATH index G1.fa\n$BWA_PATH index G2.fa\n";
my $G1_depth_sh = "$SAMTOOLS_PATH depth -aa";
my $G2_depth_sh = "$SAMTOOLS_PATH depth -aa";
$sit = 0;
for (;$sit < $G1_n;$sit++){
print "$BWA_PATH mem -t $THREAD G1.fa $G1_1[$sit] $G1_2[$sit] | $SAMTOOLS_PATH view -bS - > ./G1_BWA/G1_$sit.bam\n";
print "$SAMTOOLS_PATH sort -o ./G1_BWA/G1_$sit.sort.bam ./G1_BWA/G1_$sit.bam\n";
print "$SAMTOOLS_PATH index ./G1_BWA/G1_$sit.sort.bam\n";
$G1_depth_sh .= " ./G1_BWA/G1_$sit.sort.bam";
}

$sit = 0;
for (;$sit < $G2_n;$sit++){
print "$BWA_PATH mem -t $THREAD G1.fa $G2_1[$sit] $G2_2[$sit] | $SAMTOOLS_PATH view -bS - > ./G1_BWA/G2_$sit.bam\n";
print "$SAMTOOLS_PATH sort -o ./G1_BWA/G2_$sit.sort.bam ./G1_BWA/G2_$sit.bam\n";
print "$SAMTOOLS_PATH index ./G1_BWA/G2_$sit.sort.bam\n";
$G1_depth_sh .= " ./G1_BWA/G2_$sit.sort.bam";
}

print "$G1_depth_sh > G1.depth\n";

$sit = 0;
for (;$sit < $G1_n;$sit++){
print "$BWA_PATH mem -t $THREAD G2.fa $G1_1[$sit] $G1_2[$sit] | $SAMTOOLS_PATH view -bS - > ./G2_BWA/G1_$sit.bam\n";
print "$SAMTOOLS_PATH sort -o ./G2_BWA/G1_$sit.sort.bam ./G2_BWA/G1_$sit.bam\n";
print "$SAMTOOLS_PATH index ./G2_BWA/G1_$sit.sort.bam\n";
$G2_depth_sh .= " ./G2_BWA/G1_$sit.sort.bam";
}

$sit = 0;
for (;$sit < $G2_n;$sit++){
print "$BWA_PATH mem -t $THREAD G2.fa $G2_1[$sit] $G2_2[$sit] | $SAMTOOLS_PATH view -bS - > ./G2_BWA/G2_$sit.bam\n";
print "$SAMTOOLS_PATH sort -o ./G2_BWA/G2_$sit.sort.bam ./G2_BWA/G2_$sit.bam\n";
print "$SAMTOOLS_PATH index ./G2_BWA/G2_$sit.sort.bam\n";
$G2_depth_sh .= " ./G2_BWA/G2_$sit.sort.bam";
}

print "$G2_depth_sh > G2.depth\n";

print "perl $INSTALL_PATH" . "/step2.pl G1.fa G2.fa $G1_n $G2_n G1.depth G2.depth $CONTIG_SPECIFIC_RATIO $GAP_LENGTH $MIX_LENGTH\n";

