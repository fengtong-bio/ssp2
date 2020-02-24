#!/usr/bin/perl -w 
#
#Author: Feng Tong
#
open G1FA,"$ARGV[0]" || die;
open G2FA,"$ARGV[1]" || die;
open G1DEPTH,"$ARGV[4]" || die;
open G2DEPTH,"$ARGV[5]" || die;
my $G1_number = $ARGV[2];
my $G2_number = $ARGV[3];
my $contig_specific_ratio = $ARGV[6];
my $gap_length = $ARGV[7];
my $mix_length = $ARGV[8];
open G1DE_OUT,">G1_depth.out" || die;
open G2DE_OUT,">G2_depth.out" || die;
open G1LE_OUT,">G1_depth_length.out" || die;
open G2LE_OUT,">G2_depth_length.out" || die;
open G1FA_OUT,">G1_specific.fa" || die;
open G2FA_OUT,">G2_specific.fa" || die;

my (%G1CTG,%G2CTG,$now,@depth,$depth_sit,$real_ratio,$depth_sit2);
while (<G1DEPTH>){
	chomp;
	$now = $_;
	@depth = split(/\t/,$now);
	$depth_sit = 2;
	my $G1_count = 0;
	my $G2_count = 0;
	for (;$depth_sit <= $#depth;$depth_sit++){
		$depth_sit2 = $depth_sit - 1;
		if ($depth_sit2 <= $G1_number){
			if ($depth[$depth_sit] > 0){
				$G1_count++;
			}
		}else{
			if ($depth[$depth_sit] > 0){
				$G2_count++;
			}
		}
	}
	$real_ratio = $G1_count/$G1_number;
	if ($real_ratio >= $contig_specific_ratio && $G2_count == 0){
		print G1DE_OUT "$now\n";
		$G1CTG{$depth[0]} .= "-$depth[1]";
	}
}

while (<G2DEPTH>){
	chomp;
	$now = $_;
	@depth = split(/\t/,$now);
	$depth_sit = 2;
	my $G1_count = 0;
	my $G2_count = 0;
	for (;$depth_sit <= $#depth;$depth_sit++){
		$depth_sit2 = $depth_sit - 1;
		if ($depth_sit2 <= $G1_number){
			if ($depth[$depth_sit] > 0){
				$G1_count++;
			}
		}else{
			if ($depth[$depth_sit] > 0){
				$G2_count++;
			}
		}
	}
	$real_ratio = $G2_count/$G2_number;
	if ($real_ratio >= $contig_specific_ratio && $G1_count == 0){
		print G2DE_OUT "$now\n";
		$G2CTG{$depth[0]} .= "-$depth[1]";
	}
}

my ($start,$last,$last_gap,$long,$end,%G1L,%G2L);
foreach (keys %G1CTG){
	$now = $_;
	my @g1c = split(/-/,$G1CTG{$now});
	my $g1c_sit = 1;
	for (;$g1c_sit <= $#g1c;$g1c_sit++){
		$end = $g1c[$g1c_sit];
		if ($g1c_sit == 1){
			$start = $g1c[$g1c_sit];
			$last = $g1c[$g1c_sit];
			next;
		}
		$last_gap = $g1c[$g1c_sit] - $last;
		if ($last_gap <= $gap_length){
			$last = $g1c[$g1c_sit];
		}else{
			$long = $last - $start;
			print G1LE_OUT "$now\t$start\t$last\t$long\n";
			if ($long >= $mix_length){
				$G1L{$now} = "$start\t$last\t$long";
			}
			$start = $g1c[$g1c_sit];
			$last = $g1c[$g1c_sit];
		}
	}
	if ($end == $last){
		$long = $last - $start;
		print G1LE_OUT "$now\t$start\t$last\t$long\n";
		if ($long >= $mix_length){
			$G1L{$now} = "$start\t$last\t$long";
		}
	}
}

foreach (keys %G2CTG){
	$now = $_;
	my @g2c = split(/-/,$G2CTG{$now});
	my $g2c_sit = 1;
	for (;$g2c_sit <= $#g2c;$g2c_sit++){
		$end = $g2c[$g2c_sit];
		if ($g2c_sit == 1){
			$start = $g2c[$g2c_sit];
			$last = $g2c[$g2c_sit];
			next;
		}
		$last_gap = $g2c[$g2c_sit] - $last;
		if ($last_gap <= $gap_length){
			$last = $g2c[$g2c_sit];
		}else{
			$long = $last - $start;
			print G2LE_OUT "$now\t$start\t$last\t$long\n";
			if ($long >= $mix_length){
				$G2L{$now} = "$start\t$last\t$long";
			}
			$start = $g2c[$g2c_sit];
			$last = $g2c[$g2c_sit];
		}
	}
	if ($end == $last){
		$long = $last - $start;
		print G2LE_OUT "$now\t$start\t$last\t$long\n";
		if ($long >= $mix_length){
			$G2L{$now} = "$start\t$last\t$long";
		}
	}
}

my $control = 0;
my ($now_get,%G1F,%G2F,@G,$name,$seq);
while (<G1FA>){
chomp;
$now = $_;
	if ($now =~ />/){
		$now =~ />(\S+)/;
		$now_get = $1;
		if (exists $G1L{$now_get}){
			$control = 1;
			$name = $now_get;
		}else{
			$control = 0;
			$name = ();
		}
	}else{
		if ($control == 1){
			$G1F{$name} .= "$now";
		}
	}
}

foreach (keys %G1F){
	$now = $_;
	@G = split(/\t/,$G1L{$now});
	$seq = substr($G1F{$now},$G[0],$G[2]);
	print G1FA_OUT ">$now&$G[0]-$G[1]\n$seq\n";
}

while (<G2FA>){
	chomp;
	$now = $_;
	if ($now =~ />/){
		$now =~ />(\S+)/;
		$now_get = $1;
		if (exists $G2L{$now_get}){
			$control = 1;
			$name = $now_get;
		}else{
			$control = 0;
			$name = ();
		}
	}else{
		if ($control == 1){
			$G2F{$name} .= "$now";
		}
	}
}

foreach (keys %G2F){
	$now = $_;
	@G = split(/\t/,$G2L{$now});
	$seq = substr($G2F{$now},$G[0],$G[2]);
	print G2FA_OUT ">$now&$G[0]-$G[1]\n$seq\n";
}


