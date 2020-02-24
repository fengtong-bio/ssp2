#!/usr/bin/perl -w 
#
#Author: Feng Tong
#
my ($now,$now_path,$work_path,@G1_1,@G1_2,@G2_1,@G2_2,$specific_kmer_ratio);
my $G1_n = 0;
my $G2_n = 0;
open CFG,"$ARGV[0]" || die;
while (<CFG>){
	chomp;
	$now = $_;
	if ($now =~ /G1_1_PATH=/){
		$now =~ /G1_1_PATH=(\S+)/;
		$now_path = $1;
		$G1_1[$G1_n] = $now_path;
	}elsif ($now =~ /G1_2_PATH=/){
		$now =~ /G1_2_PATH=(\S+)/;
		$now_path = $1;
		$G1_2[$G1_n] = $now_path;
		$G1_n++;
	}elsif ($now =~ /G2_1_PATH=/){
		$now =~ /G2_1_PATH=(\S+)/;
		$now_path = $1;
		$G2_1[$G2_n] = $now_path;
	}elsif ($now =~ /G2_2_PATH=/){
		$now =~ /G2_2_PATH=(\S+)/;
		$now_path = $1;
		$G2_2[$G2_n] = $now_path;
		$G2_n++;
	}elsif ($now =~ /SPECIFIC_KMER_RATIO=/){
		$now =~ /SPECIFIC_KMER_RATIO=(\S+)/;
		$specific_kmer_ratio = $1;
	}elsif ($now =~ /WORK_PATH=/){
		$now =~ /WORK_PATH=(\S+)/;
		$work_path = $1;
	}
}
my $G1_n2 = $G1_n * 2;
my (@a,$sit,$sit_1,$G1_specific,$G2_specific,%G1,%G2);
open ALL,"$ARGV[1]" || die;
open G1_K,">$work_path/G1_specific.out" || die;
open G2_K,">$work_path/G2_specific.out" || die;
while (<ALL>){
	chomp;
	@a = split;
	my ($kmer_number,$kmer_name);
	my $G1_nn = 0;
	my $G2_nn = 0;
	$sit = 0;
	for (;$sit <= $#a;$sit += 2){
		if ($a[$sit] =~ /-/){
			$kmer_number .= "\t0";
			if ($sit < $G1_n2){
				$G1_nn++;
			}else{
				$G2_nn++;
			}
		}else{
			$sit_1 = $sit + 1;
			$kmer_number .= "\t$a[$sit_1]";
			$kmer_name = $a[$sit];
		}
	}

	my $G1_ratio = $G1_nn/$G1_n;
	my $G2_ratio = $G2_nn/$G2_n;
	if ($G1_ratio >= $specific_kmer_ratio && $G2_nn == 0){
		print G2_K "$kmer_name$kmer_number\n";
		$G2{$kmer_name} = 0;
	}elsif ($G2_ratio >= $specific_kmer_ratio && $G1_nn == 0){
		print G1_K "$kmer_name$kmer_number\n";
		$G1{$kmer_name} = 0;
	}
}


my $extract = 0;
my (%G1_paired,$tr_now,$now_n,$now_sit,$seq);
my $G1_paired_n = 1;
my $G1_paired_sit = 1;
for (;$extract < $G1_n;$extract++){
	open G1R1,"gzip -dc $G1_1[$extract]|" || die;
	open G1R2,"gzip -dc $G1_2[$extract]|" || die;
	while (<G1R1>){
		chomp;
		$now = $_;
		if ($G1_paired_sit == 1){
			$G1_paired_sit++;
			next;
		}
		if ($G1_paired_sit == 2){
			$G1_paired_sit++;
			$tr_now = reverse $now;
			$tr_now =~ tr/atcguATCGU/tagcaTAGCA/;
			$now_sit = 0;
			$now_n = length $now;
			$now_n = $now_n - 20;
			for (;$now_sit < $now_n;$now_sit++){
				$seq = substr($now,$now_sit,21);
				if (exists $G1{$seq}){
					$G1_paired{$G1_paired_n} = 0;
				}
			}
		
			$now_sit = 0;
			for (;$now_sit < $now_n;$now_sit++){
				$seq = substr($tr_now,$now_sit,21);
				if (exists $G1{$seq}){
					$G1_paired{$G1_paired_n} = 0;
				}
			}
	
			$G1_paired_n++;
			next;
		}
		if ($G1_paired_sit == 3){
			$G1_paired_sit++;
			next;
		}
		if ($G1_paired_sit == 4){
			$G1_paired_sit = 1;
			next;
		}
	}
	
	$G1_paired_n = 1;
	while (<G1R2>){
		chomp;
		$now = $_;
		if ($G1_paired_sit == 1){
			$G1_paired_sit++;
			next;
		}
		if ($G1_paired_sit == 2){
			$G1_paired_sit++;
			$tr_now = reverse $now;
			$tr_now =~ tr/atcguATCGU/tagcaTAGCA/;
			
			$now_sit = 0;
			$now_n = length $now;
			$now_n = $now_n - 20;
			for (;$now_sit < $now_n;$now_sit++){
				$seq = substr($now,$now_sit,21);
				if (exists $G1{$seq}){
					$G1_paired{$G1_paired_n} = 0;
				}
			}
			
				$now_sit = 0;
			for (;$now_sit < $now_n;$now_sit++){
				$seq = substr($tr_now,$now_sit,21);
				if (exists $G1{$seq}){
					$G1_paired{$G1_paired_n} = 0;
				}
			}
		
			$G1_paired_n++;
			next;
		}
		if ($G1_paired_sit == 3){
			$G1_paired_sit++;
			next;
		}
		if ($G1_paired_sit == 4){
			$G1_paired_sit = 1;
			next;
		}
	}
	close G1R1;
	close G1R2;
	open G1R1,"gzip -dc $G1_1[$extract]|" || die;
	open G1R2,"gzip -dc $G1_2[$extract]|" || die;
	open G1R1O,"| gzip >$work_path/kmer_reads/G1_$extract.R1.fq.gz" || die;
	open G1R2O,"| gzip >$work_path/kmer_reads/G1_$extract.R2.fq.gz" || die;
	$G1_paired_n = 1;
	while (<G1R1>){
		chomp;
		$now = $_;
		if ($G1_paired_sit == 1){
			$G1_paired_sit++;
			if (exists $G1_paired{$G1_paired_n}){
				print G1R1O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 2){
			$G1_paired_sit++;
			if (exists $G1_paired{$G1_paired_n}){
				print G1R1O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 3){
			$G1_paired_sit++;
			if (exists $G1_paired{$G1_paired_n}){
				print G1R1O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 4){
			$G1_paired_sit = 1;
			if (exists $G1_paired{$G1_paired_n}){
				print G1R1O "$now\n";
			}
			$G1_paired_n++;
		next;
		}
	}
	
	$G1_paired_n = 1;
	while (<G1R2>){
		chomp;
		$now = $_;
		if ($G1_paired_sit == 1){
			$G1_paired_sit++;
			if (exists $G1_paired{$G1_paired_n}){
				print G1R2O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 2){
			$G1_paired_sit++;
			if (exists $G1_paired{$G1_paired_n}){
				print G1R2O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 3){
			$G1_paired_sit++;
			if (exists $G1_paired{$G1_paired_n}){
				print G1R2O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 4){
			$G1_paired_sit = 1;
			if (exists $G1_paired{$G1_paired_n}){
				print G1R2O "$now\n";
			}
			$G1_paired_n++;
		next;
		}
	}
	close G1R1;
	close G1R2;
	close G1R1O;
	close G1R2O;
	
}


$extract = 0;
my %G2_paired;
$G1_paired_n = 1;
$G1_paired_sit = 1;
for (;$extract < $G2_n;$extract++){
	open G2R1,"gzip -dc $G2_1[$extract]|" || die;
	open G2R2,"gzip -dc $G2_2[$extract]|" || die;
	while (<G2R1>){
		chomp;
		$now = $_;
		if ($G1_paired_sit == 1){
			$G1_paired_sit++;
			next;
		}
		if ($G1_paired_sit == 2){
			$G1_paired_sit++;
			$tr_now = reverse $now;
			$tr_now =~ tr/atcguATCGU/tagcaTAGCA/;
			
			$now_sit = 0;
			$now_n = length $now;
			$now_n = $now_n - 20;
			for (;$now_sit < $now_n;$now_sit++){
				$seq = substr($now,$now_sit,21);
				if (exists $G2{$seq}){
					$G2_paired{$G1_paired_n} = 0;
				}
			}
		
			$now_sit = 0;
			for (;$now_sit < $now_n;$now_sit++){
				$seq = substr($tr_now,$now_sit,21);
				if (exists $G2{$seq}){
					$G2_paired{$G1_paired_n} = 0;
				}
			}
			
			$G1_paired_n++;
			next;
		}
		if ($G1_paired_sit == 3){
			$G1_paired_sit++;
			next;
		}
		if ($G1_paired_sit == 4){
			$G1_paired_sit = 1;
			next;
		}
	}
	
	$G1_paired_n = 1;
	while (<G2R2>){
		chomp;
		$now = $_;
		if ($G1_paired_sit == 1){
			$G1_paired_sit++;
			next;
		}
		if ($G1_paired_sit == 2){
			$G1_paired_sit++;
			$tr_now = reverse $now;
			$tr_now =~ tr/atcguATCGU/tagcaTAGCA/;
			
			$now_sit = 0;
			$now_n = length $now;
			$now_n = $now_n - 20;
			for (;$now_sit < $now_n;$now_sit++){
				$seq = substr($now,$now_sit,21);
				if (exists $G2{$seq}){
					$G2_paired{$G1_paired_n} = 0;
				}
			}
			
			$now_sit = 0;
			for (;$now_sit < $now_n;$now_sit++){
				$seq = substr($tr_now,$now_sit,21);
				if (exists $G2{$seq}){
					$G2_paired{$G1_paired_n} = 0;
				}
			}
		
			$G1_paired_n++;
			next;
		}
		if ($G1_paired_sit == 3){
			$G1_paired_sit++;
			next;
		}
		if ($G1_paired_sit == 4){
			$G1_paired_sit = 1;
			next;
		}
	}
	close G2R1;
	close G2R2;
	
	open G2R1,"gzip -dc $G2_1[$extract]|" || die;
	open G2R2,"gzip -dc $G2_2[$extract]|" || die;
	open G2R1O,"| gzip >$work_path/kmer_reads/G2_$extract.R1.fq.gz" || die;
	open G2R2O,"| gzip >$work_path/kmer_reads/G2_$extract.R2.fq.gz" || die;
	$G1_paired_n = 1;
	while (<G2R1>){
		chomp;
		$now = $_;
		if ($G1_paired_sit == 1){
			$G1_paired_sit++;
			if (exists $G2_paired{$G1_paired_n}){
				print G2R1O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 2){
			$G1_paired_sit++;
			if (exists $G2_paired{$G1_paired_n}){
				print G2R1O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 3){
			$G1_paired_sit++;
			if (exists $G2_paired{$G1_paired_n}){
				print G2R1O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 4){
			$G1_paired_sit = 1;
			if (exists $G2_paired{$G1_paired_n}){
				print G2R1O "$now\n";
			}
			$G1_paired_n++;
		next;
		}
	}
	
	$G1_paired_n = 1;
	while (<G2R2>){
		chomp;
		$now = $_;
		if ($G1_paired_sit == 1){
			$G1_paired_sit++;
			if (exists $G2_paired{$G1_paired_n}){
				print G2R2O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 2){
			$G1_paired_sit++;
			if (exists $G2_paired{$G1_paired_n}){
				print G2R2O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 3){
			$G1_paired_sit++;
			if (exists $G2_paired{$G1_paired_n}){
				print G2R2O "$now\n";
			}
		next;
		}
		if ($G1_paired_sit == 4){
			$G1_paired_sit = 1;
			if (exists $G2_paired{$G1_paired_n}){
				print G2R2O "$now\n";
			}
			$G1_paired_n++;
		next;
		}
	}
	close G2R1;
	close G2R2;
	close G2R1O;
	close G2R2O;
	
}
