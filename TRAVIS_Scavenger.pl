#!/usr/bin/env perl
use warnings;
use strict;
use File::Path qw(make_path);	
use File::Basename;
use File::Spec;
use Data::Dumper;
use List::MoreUtils qw(uniq); # maybe replace?

unless (@ARGV){
	die "usage: perl TRAVIS_Scavenger.pl <TRAVIS_control_center.csv>\n";
}
my $TCC_file = shift @ARGV;


print '#' x 70,"\n";
print '#' x 70,"\n";
print "\t\t\tThis is TRAVIS Scavenger v20190209\n";
print '#' x 70,"\n";
print '#' x 70,"\n";
print "\n";

#~ my $TCC_file = '/home/skaefer/sync/work/reoviridae/reo_full/reo_full_TCC.csv';
#~ my $TCC_file = '/home/skaefer/sync/work/reoviridae/reo_mini_test/reo_mini_test_TCC.csv';
my $partially ; #give a value if you want to run only the first entries of each file


###from http://lsrtools.1apps.com/wavetorgb/index.asp 
#hexcodes that represent wavelength [nm]
my %rainbow  = ( 
	401 => '#4900AA', 	402 => '#4900AE', 	403 => '#4900B3', 	404 => '#4900B7', 	405 => '#4900BC', 	406 => '#4800C0', 	407 => '#4800C4', 	408 => '#4700C9', 	409 => '#4600CD', 	410 => '#4600D2', 	411 => '#4500D6', 	412 => '#4400DB', 	413 => '#4300DF', 	414 => '#4100E4', 	415 => '#4000E8', 	416 => '#3F00ED', 	417 => '#3D00F1', 	418 => '#3C00F6', 	419 => '#3A00FA', 	420 => '#3800FF', 	421 => '#3500FF', 	422 => '#3300FF', 	423 => '#3000FF', 	424 => '#2D00FF', 	425 => '#2A00FF', 	426 => '#2700FF', 	427 => '#2400FF', 	428 => '#2200FF', 	429 => '#1F00FF', 	430 => '#1C00FF', 	431 => '#1900FF', 	432 => '#1600FF', 	433 => '#1300FF', 	434 => '#1100FF', 	435 => '#0E00FF', 	436 => '#0B00FF', 	437 => '#0800FF', 	438 => '#0500FF', 	439 => '#0200FF', 	440 => '#0000FF', 	441 => '#0005FF', 	442 => '#000AFF', 	443 => '#000FFF', 	444 => '#0014FF', 	445 => '#0019FF', 	446 => '#001EFF', 	447 => '#0023FF', 	448 => '#0028FF', 	449 => '#002DFF',
	450 => '#0033FF', 	451 => '#0038FF', 	452 => '#003DFF', 	453 => '#0042FF', 	454 => '#0047FF', 	455 => '#004CFF', 	456 => '#0051FF', 	457 => '#0056FF', 	458 => '#005BFF', 	459 => '#0060FF', 	460 => '#0066FF', 	461 => '#006BFF', 	462 => '#0070FF', 	463 => '#0075FF', 	464 => '#007AFF', 	465 => '#007FFF', 	466 => '#0084FF', 	467 => '#0089FF', 	468 => '#008EFF', 	469 => '#0093FF', 	470 => '#0099FF', 	471 => '#009EFF', 	472 => '#00A3FF', 	473 => '#00A8FF', 	474 => '#00ADFF', 	475 => '#00B2FF', 	476 => '#00B7FF', 	477 => '#00BCFF', 	478 => '#00C1FF', 	479 => '#00C6FF', 	480 => '#00CCFF', 	481 => '#00D1FF', 	482 => '#00D6FF', 	483 => '#00DBFF', 	484 => '#00E0FF', 	485 => '#00E5FF', 	486 => '#00EAFF', 	487 => '#00EFFF', 	488 => '#00F4FF', 	489 => '#00F9FF', 	490 => '#00FFFF', 	491 => '#00FFF2', 	492 => '#00FFE5', 	493 => '#00FFD8', 	494 => '#00FFCC', 	495 => '#00FFBF', 	496 => '#00FFB2', 	497 => '#00FFA5', 	498 => '#00FF99', 	499 => '#00FF8C',
	500 => '#00FF7F', 	501 => '#00FF72', 	502 => '#00FF66', 	503 => '#00FF59', 	504 => '#00FF4C', 	505 => '#00FF3F', 	506 => '#00FF33', 	507 => '#00FF26', 	508 => '#00FF19', 	509 => '#00FF0C', 	510 => '#00FF00', 	511 => '#03FF00', 	512 => '#07FF00', 	513 => '#0AFF00', 	514 => '#0EFF00', 	515 => '#12FF00', 	516 => '#15FF00', 	517 => '#19FF00', 	518 => '#1DFF00', 	519 => '#20FF00', 	520 => '#24FF00', 	521 => '#28FF00', 	522 => '#2BFF00', 	523 => '#2FFF00', 	524 => '#33FF00', 	525 => '#36FF00', 	526 => '#3AFF00', 	527 => '#3DFF00', 	528 => '#41FF00', 	529 => '#45FF00', 	530 => '#48FF00', 	531 => '#4CFF00', 	532 => '#50FF00', 	533 => '#53FF00', 	534 => '#57FF00', 	535 => '#5BFF00', 	536 => '#5EFF00', 	537 => '#62FF00', 	538 => '#66FF00', 	539 => '#69FF00', 	540 => '#6DFF00', 	541 => '#70FF00', 	542 => '#74FF00', 	543 => '#78FF00', 	544 => '#7BFF00', 	545 => '#7FFF00', 	546 => '#83FF00', 	547 => '#86FF00', 	548 => '#8AFF00', 	549 => '#8EFF00',
	550 => '#91FF00', 	551 => '#95FF00', 	552 => '#99FF00', 	553 => '#9CFF00', 	554 => '#A0FF00', 	555 => '#A3FF00', 	556 => '#A7FF00', 	557 => '#ABFF00', 	558 => '#AEFF00', 	559 => '#B2FF00', 	560 => '#B6FF00', 	561 => '#B9FF00', 	562 => '#BDFF00', 	563 => '#C1FF00', 	564 => '#C4FF00', 	565 => '#C8FF00', 	566 => '#CCFF00', 	567 => '#CFFF00', 	568 => '#D3FF00', 	569 => '#D6FF00', 	570 => '#DAFF00', 	571 => '#DEFF00', 	572 => '#E1FF00', 	573 => '#E5FF00', 	574 => '#E9FF00', 	575 => '#ECFF00', 	576 => '#F0FF00', 	577 => '#F4FF00', 	578 => '#F7FF00', 	579 => '#FBFF00', 	580 => '#FFFF00', 	581 => '#FFFB00', 	582 => '#FFF700', 	583 => '#FFF300', 	584 => '#FFEF00', 	585 => '#FFEB00', 	586 => '#FFE700', 	587 => '#FFE300', 	588 => '#FFDF00', 	589 => '#FFDB00', 	590 => '#FFD700', 	591 => '#FFD300', 	592 => '#FFCF00', 	593 => '#FFCC00', 	594 => '#FFC800', 	595 => '#FFC400', 	596 => '#FFC000', 	597 => '#FFBC00', 	598 => '#FFB800', 	599 => '#FFB400',
	600 => '#FFB000', 	601 => '#FFAC00', 	602 => '#FFA800', 	603 => '#FFA400', 	604 => '#FFA000', 	605 => '#FF9C00', 	606 => '#FF9900', 	607 => '#FF9500', 	608 => '#FF9100', 	609 => '#FF8D00', 	610 => '#FF8900', 	611 => '#FF8500', 	612 => '#FF8100', 	613 => '#FF7D00', 	614 => '#FF7900', 	615 => '#FF7500', 	616 => '#FF7100', 	617 => '#FF6D00', 	618 => '#FF6900', 	619 => '#FF6600', 	620 => '#FF6200', 	621 => '#FF5E00', 	622 => '#FF5A00', 	623 => '#FF5600', 	624 => '#FF5200', 	625 => '#FF4E00', 	626 => '#FF4A00', 	627 => '#FF4600', 	628 => '#FF4200', 	629 => '#FF3E00', 	630 => '#FF3A00', 	631 => '#FF3600', 	632 => '#FF3300', 	633 => '#FF2F00', 	634 => '#FF2B00', 	635 => '#FF2700', 	636 => '#FF2300', 	637 => '#FF1F00', 	638 => '#FF1B00', 	639 => '#FF1700', 	640 => '#FF1300', 	641 => '#FF0F00', 	642 => '#FF0B00', 	643 => '#FF0700', 	644 => '#FF0300', 	645 => '#FF0000'
);

my $color_resolution = 5; #rainbow scale will be applied to the sample ORF in chunks of this percentage
my $color_sensitivity= 0.00001; #is the blast-evalue threshold for the matches

print '#' x 70,"\n";
print '#' x 70,"\n";
print "\t\t\tThis is TRAVIS Scavenger v0.8a\n";
print '#' x 70,"\n";
print '#' x 70,"\n";
my @supported_search_tools = ('blastp', 'hmmer', 'jackhmmer', 'mmseqs', 'blastp_vs_nr'); #only entries with the supportet search tools will be evaluated

###reading TRAVIS Control Center
my %TCC; #contains configuration parameters from $TCC_file 
&read_TCC(
	\$TCC_file, #file path
	\%TCC #undef hash
);
print "connected to TRAVIS Control Center\n";
&check_dir(\$TCC{'reference_gbx'}); #genebank xmls will be stored here locally


my $sample_library =$TCC{'sample_library'};


my %sample_library_entries;
my @samples;
my @gbx_times;
open (my $SAMPLE_LIBRARY, '<', $sample_library) or die "could not read from '$sample_library': $!\n";
while (my $line = <$SAMPLE_LIBRARY>){
	chomp $line;
	#~ print $line, "\n";
	my @cols = split(',', $line);
	unless ($sample_library_entries{'header'}){
		#~ print "reading header\n";
		@{$sample_library_entries{'header'}} = @cols;
	} else {
		
		#~ print "saving data for $cols[1]\n";
		push (@samples, $cols[1]);
		@{$sample_library_entries{$cols[1]}} = @cols
	}
}
push (@{$sample_library_entries{'header'}} ,'suspicious_sequences');




#already evaluated sample IDs will be stored here so that you can resume evaluation
#resume calculation can be switched on/off in TCC
my $scavenged_samples = File::Spec -> catfile ($TCC{'result_dir'}, 'scavenged_samples.txt');
my %already_scavenged;
if ((-s $scavenged_samples) and ($TCC{'resume_calculation'})){
	print "will try to resume calculation!\n";
	open (my $SCAVENGED_SAMPLES, '<', $scavenged_samples) or die "could not read from '$scavenged_samples': $!\n";
	while (my $line = <$SCAVENGED_SAMPLES>){
		chomp $line;
		print "already scavenged $line\n";
		$already_scavenged{$line} = 1;
	}
	close($SCAVENGED_SAMPLES);
}
open (my $SCAVENGED_SAMPLES, '>', $scavenged_samples) or die "could not write to '$scavenged_samples': $!\n";


my $all_NT_col; #will contain the index of the column called 'All_NT_ACC'
my $main_NT_col; #will contain the index of the column called '<MAIN>_NT_ACC'
my $main_PID_col; #will contain the index of the column called '<MAIN>_PID'
my @main_PIDs;
my $main_name; #will contain thename of the main protein/gene/whatever '<MAIN>_NT_ACC'


#reading the preparsed reference library
(my $preparsed_reference_library = $TCC{'reference_library'}) =~ s/\.csv$/_preparsed.csv/;



if (-s $preparsed_reference_library){
	open (my $REFLIB, '<', $preparsed_reference_library) or die "could not read from '$preparsed_reference_library': $!\n";
	print "reading reference library '$preparsed_reference_library'...\n";
	while (my $line = <$REFLIB>){
		chomp $line;
		my @cols = split (',', $line);
		###check reference library table headers
		unless ($main_PID_col){ #do this if $main_PID_col is not yet defined
			###get the accession number column numbers by looping through them
			#example: All_NT_ACC,RdRp_NT_ACC,RdRp_PID
			for (my $i=0; $i < scalar @cols; $i++){
				if ($cols[$i] eq 'All_NT_ACC'){
					$all_NT_col= $i;
				} elsif ($cols[$i] =~ m/(.*)_NT_ACC/){
					$main_NT_col= $i;
					$main_name = $1;
				} elsif (($main_NT_col) and ($cols[$i] =~ m/$main_name\_PID/)){
					$main_PID_col= $i;
				}
			}
		} else {
			(my $crt_PID = $cols[$main_PID_col])  =~ s/\.\d\s*$//; #get rid of version number
			push (@main_PIDs, $crt_PID);
		}
	}
} else {
	print "this seems to be a blind flight...let's see....\n";
	$TCC{'reference_library'} = 'blind' unless $TCC{'reference_library'};
}



my %suspicious_sequences; #1st dimension: Sample ID, 2nd dimension: sequence header, 3rd dimension: ORF#, 4th dimension: methods that identified the sequence, 5th dimension: possible annotation
my %infection_status; #infected/clean of Sample ID
my $ORF_results = File::Spec -> catfile ($TCC{'result_dir'},  'ORF_results.csv');
print "reading TRAVIS Core results '$ORF_results'...\n";
&get_times(\$ORF_results);
open (my $ORF_RESULTS, '<', $ORF_results) or die "could not read from '$ORF_results': $!\n";
print '#' x 70,"\n";
my @infected; #will collect all sample IDs of infected samples
###reading the results and sorting according to sample, search tool/method and the detected related sequence(s)
my %reflib_matches_counter;
my %nr_matches_counter;
my $entrycounter;
while (my $line = <$ORF_RESULTS>){
	chomp $line;
	#~ next if ($line =~ m/^#/);
	
	# last if ($entrycounter >= 30);
	my @cols = split (',', $line);
	if ($cols[1] ne 'NA'){
		$cols[1] =~ m/(.+)_(ORF_\d+)/;
		
		my $seq_pt1 = $1;
		my $seq_pt2 = $2;
		
		if ($cols[3] eq 'nr'){
			$entrycounter++;
			$nr_matches_counter{$cols[1]}++;
			next if ($nr_matches_counter{$cols[1]} > $TCC{'max_references'} );
		} else {
			next if ($TCC{'reference_library'} eq 'blind');
			$reflib_matches_counter{$cols[1]}{$cols[3]}++;
			next if ($reflib_matches_counter{$cols[1]}{$cols[3]} > $TCC{'max_references'} );
		}
		
		
		if ($cols[2] =~ m/acc\|/){
			my @acc_entries = split ('\|', $cols[2]);
			die "fuck $cols[2] if\n" if !(@acc_entries);
			$cols[2] = $acc_entries[2];
			$cols[2] =~ s/\.\d$//;
			$cols[2] .= '__NA';
			print "extracted accession number $cols[2]\n";
		
		}
		
		
		push ( @{$suspicious_sequences{$cols[0]}{$seq_pt1}{$seq_pt2}{$cols[5]}},$cols[3].'___'.$cols[2]);
		$infection_status{$cols[0]} = 'positive';
	} else {
		$infection_status{$cols[0]} = 'clean' unless ($infection_status{$cols[0]});
	}
	if ($infection_status{$cols[0]} eq 'positive'){
		push (@infected, $cols[0]) unless ($cols[0] ~~@infected);
	}
}

my $infected_counter;

if (@infected){
	my %main_hit_sequences; #collecting the sequences of the main set for later alignment and phylogeny reconstruction
	print scalar @infected, ' of ', scalar keys %infection_status, " samples are possibly infected:\n";
	###loop through the detailed results of all infected samples
	foreach my $sample_ID (@infected){
		if ($already_scavenged{$sample_ID}){
			print "already scavenged $sample_ID!\n";
			print $SCAVENGED_SAMPLES $sample_ID, "\n";  #remember that it has been evaluated
			next;
		}
		my $t1  = time();
		$infected_counter++;
		#~ last if ($infected_counter >=2); 
		my %ORF_related_to; #collects all related PIDs for the suspicious ORFs
		my %ORF_relatedness_to_PID_identified_by_method; #1st dim: PID, 2nd dim: ORF ID, value: array listing methods
		my %suspicious_ORFs_AA; #collects the AA sequences of the ORFs of the sample
		my %reference_ORFs_AA; #collects the AA sequences of the ORFs of the related references
		my %ORF_needs_a_rainbow; #stores combinations of sample and reference ORFs that will get a coloration
		my %found_by; #stores the information which search-tool/method identified the sample ORF
		my %start_end_orientation_by_PID; #stores the start, end and orientation of a PID, separated by underscore
		#~ my %reference_AA_to_sample_ORF; # reference_AA (NT_ACC_ORF-Nr) = PID
		# NT_ACC_ORF-Nr     header_ORF-Nr        hit-PID
		my %potential_annotations; #stores all potential annotations by various methods for the sample ORFs. seen in title of plot element
		my %sequences; #headbreaking multidimensional hash that contains plot information about aaaaaall sequences:
		#####$sequences{$group}{$sequence_header}{$sequence_type(NT or AA)}{$sequence_description}{$start_end_orientation[0]}
		#~ my %sequence_refs;
		my %check_ORF_related_to;
		my %check_ORF_relatedness_to_PID_identified_by_method;
		my %check_potential_annotations;
		my %check_found_by;
		
		print '#' x 70,"\n";
		my $number_of_suspicious = scalar keys %{$suspicious_sequences{$sample_ID}};
		print "$sample_ID: $number_of_suspicious sequence(s)!\n";
		push (@{$sample_library_entries{$sample_ID}}, $number_of_suspicious);
		my $crt_sample_info;
		for (my $header_position = 2; $header_position < scalar @{$sample_library_entries{'header'}}; $header_position++){
			$crt_sample_info .= ${$sample_library_entries{'header'}}[$header_position].': '. ${$sample_library_entries{$sample_ID}}[$header_position]."\n";
		}

		#~ print $crt_sample_info;
		my $crt_ORF_CSV = File::Spec -> catfile ($TCC{'ORF_dir'}, $sample_ID.'_ORF_data.csv');
		###loading all ORF data from the sample ORF file
		print "retrieving ORF data from '$crt_ORF_CSV'...\n";
		my %ORF_data;
		my %ORF_data_by_sequence;
		open (my $CRT_ORF_CSV, '<', $crt_ORF_CSV) or die "could not read from '$crt_ORF_CSV': $!\n";
		my $crt_header;
		while (my $line = <$CRT_ORF_CSV>){
			chomp $line;
			$line =~ m/^([^,]+),(.+)$/;
			my $crt_ORF_header = $1;
			my $crt_ORF_data = $2;
			$ORF_data{$crt_ORF_header} = $crt_ORF_data;
			if ($line =~ m/^>([^,]+)/){
				$crt_header = $1;
			} else {
				$ORF_data_by_sequence{$crt_header}{$crt_ORF_header} = $crt_ORF_data;
			}
		}
		
		###getting the suspicious sequences from that sample 
		
		foreach my $sequence_header (keys %{$suspicious_sequences{$sample_ID}}){
			print "\t\t$sequence_header\n";
			my $crt_full_seq = $sequence_header.','.$ORF_data{'>'.$sequence_header};
			foreach my $crt_ORF (sort keys %{$ORF_data_by_sequence{$sequence_header}}){
				my @cols = split (',', $ORF_data_by_sequence{$sequence_header}{$crt_ORF});
				$suspicious_ORFs_AA{$crt_ORF} = $cols[7];
			}
		}
		
		
		
		
		###loop through the TC results:
		foreach my $sequence_header (sort keys %{$suspicious_sequences{$sample_ID}}){
			
			
			###loop through all ORFs on that NT
			foreach my $crt_ORF (sort keys %{$ORF_data_by_sequence{$sequence_header}}){
				print "\t\t\tevaluating $crt_ORF\n";
				###get sequence original header and ORF number
				$crt_ORF =~ m/(.+)_(ORF_\d+)/;
				my $crt_ORF_number  = $2;
				###loop through suppored search tools
				if ( $suspicious_sequences{$sample_ID}{$sequence_header}{$crt_ORF_number}){
					foreach my $crt_method (sort keys %{$suspicious_sequences{$sample_ID}{$sequence_header}{$crt_ORF_number}}){
						foreach my $identified_as (@{$suspicious_sequences{$sample_ID}{$sequence_header}{$crt_ORF_number}{$crt_method}}){
							my @crt_triple_elements = split('___',$identified_as); #splits the file source info
							#~ print "!#!#!#!#!#!#!#!de triple elements based on $identified_as\n";
							
							unless ($check_found_by{$crt_ORF}{$crt_method}){
								push(@{$found_by{$crt_ORF}}, $crt_method);
								$check_found_by{$crt_ORF}{$crt_method} = 1;
							}
							if ($crt_method eq 'hmmer'){
								#~ print "\t\t\t\t\t\t\tit was hmmer! will search the cluster file for accessions...\n";
								#~ my @crt_triple_elements = split('___',$identified_as);
								my $crt_cluster_file = File::Spec -> catfile($TCC{'reference_fastas'}, $TCC{'database_name'}.'_'.$crt_triple_elements[0].'.fasta');
								unless (-s $crt_cluster_file ){
									$crt_cluster_file = File::Spec -> catfile($TCC{'reference_fastas'},$TCC{'database_name'}.'_AA_'.$crt_triple_elements[0].'.fasta');
								}
								unless (-s $crt_cluster_file ){ # in case of a blind flight
									$crt_cluster_file = File::Spec -> catfile($TCC{'reference_fastas'},$crt_triple_elements[0].'.fasta');
								}
								open (my $CRT_CLUSTER_FILE, '<', $crt_cluster_file) or die "could not read from '$crt_cluster_file': $!\n";
								while (my $line = <$CRT_CLUSTER_FILE>){
									chomp $line;
									if ($line =~ m/^>(.*)/){
										my @crt_double_elements = split('__', $1);
										#~ my @crt_double_elements = split('__', $crt_triple_elements[-1]);
										#~ print "\t\t\t\t\t\t\t\tthis one! $crt_double_elements[0]\n";
										$crt_double_elements[-1] =~ s/^_//;
										#~ print "hmmer says $crt_double_elements[-1]\n";
										#~ print Dumper @crt_double_elements;
										unless ($check_ORF_related_to{$crt_ORF}{$crt_double_elements[0]}){
											push(@{$ORF_related_to{$crt_ORF}}, $crt_double_elements[0]);
											$check_ORF_related_to{$crt_ORF}{$crt_double_elements[0]} = 1;
										}
										unless ( $check_ORF_relatedness_to_PID_identified_by_method{$crt_double_elements[0]}{$crt_ORF}{$crt_method}){
											push(@{$ORF_relatedness_to_PID_identified_by_method{$crt_double_elements[0]}{$crt_ORF}}, $crt_method);
											$check_ORF_relatedness_to_PID_identified_by_method{$crt_double_elements[0]}{$crt_ORF}{$crt_method} = 1;
										}
										unless ($check_potential_annotations{$crt_ORF}{$crt_double_elements[-1]}){
											push(@{$potential_annotations{$crt_ORF}}, $crt_double_elements[-1]);
											$check_potential_annotations{$crt_ORF}{$crt_double_elements[-1]} = 1;
										}
#~ push(@{$identified_by{$crt_double_elements[0]}}, $crt_method) unless ($crt_method ~~ @{$identified_by{$crt_double_elements[0]}});
												}
								}
								close ($CRT_CLUSTER_FILE);	
							} elsif ($crt_method eq 'blastp_vs_nr'){
								my @crt_double_elements = split('__', $crt_triple_elements[1]);
								#~ print "\t\t\t\t\t\t\t\tthis one! $crt_double_elements[0]\n";
								#~ print Dumper @crt_double_elements;
								#~ print "$crt_method says $crt_double_elements[-1]\n";
								
								my $crt_description;
								my $crt_sequence;
								&fetch_sequence(
									\$crt_double_elements[-1], #PID
									\$crt_description, #NCBI description of the sequence to be retrieved
									\$crt_sequence, #sequence to be retrieved
									\'protein' #which database
								);
								unless ($check_ORF_related_to{$crt_ORF}{$crt_double_elements[0]}){
									push(@{$ORF_related_to{$crt_ORF}}, $crt_double_elements[0]);
									$check_ORF_related_to{$crt_ORF}{$crt_double_elements[0]} = 1;
								}
								unless ( $check_ORF_relatedness_to_PID_identified_by_method{$crt_double_elements[0]}{$crt_ORF}{$crt_method}){
									push(@{$ORF_relatedness_to_PID_identified_by_method{$crt_double_elements[0]}{$crt_ORF}}, $crt_method) ;
									$check_ORF_relatedness_to_PID_identified_by_method{$crt_double_elements[0]}{$crt_ORF}{$crt_method} = 1;
								}
								unless ($check_potential_annotations{$crt_ORF}{$crt_description}){
									push(@{$potential_annotations{$crt_ORF}}, $crt_description);
									$check_potential_annotations{$crt_ORF}{$crt_description} = 1;
								}
								
								
							} else {
								#~ ###if it was not hmmer, there is only a specific sequence to retrieve
								#~ print "\t\t\t\t\t\t\tit was something where i should be able to find an accession number...\n";
								#~ my @crt_triple_elements = split('___',$identified_as);
								my @crt_double_elements = split('__', $crt_triple_elements[1]);
								#~ print "\t\t\t\t\t\t\t\tthis one! $crt_double_elements[0]\n";
								#~ print Dumper @crt_double_elements;
								#~ print "$crt_method says $crt_double_elements[-1]\n";
								unless ($check_ORF_related_to{$crt_ORF}{$crt_double_elements[0]}){
									push(@{$ORF_related_to{$crt_ORF}}, $crt_double_elements[0]);
									$check_ORF_related_to{$crt_ORF}{$crt_double_elements[0]} = 1;
								}
								unless ($check_ORF_relatedness_to_PID_identified_by_method{$crt_double_elements[0]}{$crt_ORF}{$crt_method}){
									push(@{$ORF_relatedness_to_PID_identified_by_method{$crt_double_elements[0]}{$crt_ORF}}, $crt_method);
									$check_ORF_relatedness_to_PID_identified_by_method{$crt_double_elements[0]}{$crt_ORF}{$crt_method} = 1;
								}
								unless ( $check_potential_annotations{$crt_ORF}{$crt_double_elements[-1]}){
									push(@{$potential_annotations{$crt_ORF}}, $crt_double_elements[-1]);
									$check_potential_annotations{$crt_ORF}{$crt_double_elements[-1]} = 1;
								}
								#~ push(@{$identified_by{$crt_double_elements[0]}}, $crt_method) unless ($crt_method ~~ @{$identified_by{$crt_double_elements[0]}});
								#~ push(@TRAVIS_annotation_PIDs, $crt_double_elements[0]) unless ($crt_double_elements[0] ~~ @TRAVIS_annotation_PIDs);
							
							}
						}
					}
				} else {
					#~ print "\t\t\t\t\t\tTRAVIS Core found no annotations for this ORF...\n";
				}	
			}
		}
	
		#again looping through all the ORFs. now with complete data!
		my $seqcountdings;
		foreach my $sequence_header (sort keys %{$suspicious_sequences{$sample_ID}}){
			print '#' x 70,"\n";
			print "current sequence basename: $sequence_header\n";
			$seqcountdings++;
			#~ last if ($seqcountdings >= 2);
			my $NT_length = length ($ORF_data{'>'.$sequence_header});
			#~ print " full NT length: $NT_length\n";
			my $crt_full_seq = $sequence_header.','.$ORF_data{'>'.$sequence_header};	
			#~ print " $crt_full_seq\n";
			$sequences{$sample_ID.'|||'.$sequence_header}{$sequence_header}{'NT'}{$sequence_header}{'1_'.$NT_length}= "Origin: $sample_ID\n"; #generating plot entry for base NT sequence
			###loop through all ORFs on that NT
			my %NC_ACC_by_PID; #stores the NT_ACCs of PIDs
			foreach my $crt_ORF (sort keys %{$ORF_data_by_sequence{$sequence_header}}){
				#~ print "current ORF: $crt_ORF\n";								
				
				my @cols = split (',', $ORF_data_by_sequence{$sequence_header}{$crt_ORF});
				#~ ###check orientation of the ORF
				my $crt_from;
				my $crt_to; 
				if ($cols[2] < $cols[3]){
					$crt_from = $cols[2];
					$crt_to = $cols[3];
				} else {
					$crt_from = $cols[3];
					$crt_to = $cols[2];
				}
				my $crt_orientation;
				if ($cols[0] =~ m/reverse/){
					$crt_orientation = 'reverse';
				} else {
					$crt_orientation = $cols[0];
				}
				#~ print "  from $crt_from to $crt_to, $crt_orientation\n";
				my $crt_ORF_start_end_orientation =  $crt_from.'_'.$crt_to.'_'.$crt_orientation;
				my $crt_ORF_found_by_methods;
				if ($found_by{$crt_ORF}){
					$ORF_needs_a_rainbow{$crt_ORF} = 'self';
					#~ print "\t found by:\n";
					#~ foreach my $crt_method (@{$found_by{$crt_ORF}}){
						#~ print "\t  $crt_method\n";
					#~ }
					$crt_ORF_found_by_methods = join("\n", sort (@{$found_by{$crt_ORF}}));
					foreach my $crt_PID (@{$ORF_related_to{$crt_ORF}}){
						print "\t\tcurrent PID: $crt_PID\n";
						my $methods = join(' & ', sort (@{$ORF_relatedness_to_PID_identified_by_method{$crt_PID}{$crt_ORF}}));
						#~ print "\t\t by methods: $methods\n";
						print "\t\t  retrieving NCBI data...\n";
						my %crt_PID_entries;
						my %crt_NT_entries;
						&fetch_gbx(\$crt_PID,\%crt_PID_entries,\'protein',\$TCC{'reference_gbx'}); #retrieves a genebank entry in xml format for the current PID
						if ($crt_PID_entries{'main'}{'source-db'}){ #and checks, if it finds a NT accession number as the source/origin of the PID and if there is one, it will
							my @associated_PIDs;
							#~ $crt_entries{'main'}{'source-db'} =~ m/([^\s]+)\.*\d*$/;
							$crt_PID_entries{'main'}{'source-db'} =~ m/([A-Z][A-Z]+_*\d{3}\d*)\.*\d*$/;
							my $crt_NT_ACC;
							if ($1) {
								$crt_NT_ACC = $1;
								$crt_NT_ACC =~ s/\.\d+\s*$//;
								
								&fetch_gbx(\$crt_NT_ACC,\%crt_NT_entries,\'nuccore',\$TCC{'reference_gbx'}); #retrieve the genebank entry of the source 
								if (%crt_NT_entries){
									#~ print "\t\t\thas its origin in $crt_NT_ACC\n";
									$sequences{$sample_ID.'|||'.$sequence_header}{$crt_NT_ACC}{'NT'}{$crt_NT_ACC.'_'.$crt_NT_entries{'main'}{'definition'}}{'1_'.$crt_NT_entries{'main'}{'length'}}= 'Taxonomy: '.$crt_NT_entries{'main'}{'taxonomy'};

									#~ print "got some nucleotide stuff!\n";
									#~ print "\t\t\t looking for annotated ORFs...\n";
									my %ORF_names; #orf names by start_end_orientation
									foreach my $entry (sort keys %crt_NT_entries){ #looks through the features if there are more annotated ORFs that can be used for comparison and plotting in this reference
										next if (!$crt_NT_entries{$entry}{'from'} or !$crt_NT_entries{$entry}{'to'} or !$crt_NT_entries{$entry}{'translation'});
										
										
										if ($crt_NT_entries{$entry}{'product'}){
											$crt_NT_entries{$entry}{'product'} =~ s/[^a-zA-Z0-9_\-]/_/g;
											$crt_NT_entries{$entry}{'protein_id'} =~ s/\.\d+\s*$//;
											push (@associated_PIDs ,$crt_NT_entries{$entry}{'protein_id'});
											###checks orientation of the ORF
											if ($crt_NT_entries{$entry}{'from'} > $crt_NT_entries{$entry}{'to'}){
												$start_end_orientation_by_PID{$crt_NT_entries{$entry}{'protein_id'}} = $crt_NT_entries{$entry}{'to'}.'_'.$crt_NT_entries{$entry}{'from'}.'_reverse';
											} else {
												$start_end_orientation_by_PID{$crt_NT_entries{$entry}{'protein_id'}}= $crt_NT_entries{$entry}{'from'}.'_'.$crt_NT_entries{$entry}{'to'}.'_forward';
											}
											$ORF_names{$start_end_orientation_by_PID{$crt_NT_entries{$entry}{'protein_id'}}} = $crt_NT_entries{$entry}{'protein_id'}.'   '.$crt_NT_entries{$entry}{'product'} ;
											
											my $crt_reference_ORF = $crt_NT_ACC.'_ORF_'.$crt_NT_entries{$entry}{'protein_id'};
											###if the PID of the ORF is the PID that is the match, it will be marked for generating a rainbow
											if ($crt_PID eq $crt_NT_entries{$entry}{'protein_id'}){
												#~ print "setting refernence for $sequence_header.'_rainbow_'.$crt_reference_ORF\n";
												$ORF_needs_a_rainbow{$sequence_header.'_rainbow_'.$crt_reference_ORF} = $crt_ORF;
												#~ $ORF_needs_a_rainbow{$crt_reference_ORF} = $crt_ORF;
												$reference_ORFs_AA{$crt_reference_ORF} = $crt_NT_entries{$entry}{'translation'};
												#~ print "\t\t\t  running rainbow on $crt_reference_ORF vs $crt_ORF\n" ;
												#~ print "\t\t\t\tthat would be:\n\t\t\t\t$crt_ORF: $suspicious_ORFs_AA{$crt_ORF}\n";
												#~ print "\t\t\t\tvs.\n\t\t\t\t$crt_reference_ORF: $reference_ORFs_AA{$crt_reference_ORF}\n";
												$sequences{$sample_ID.'|||'.$sequence_header}{$crt_NT_ACC}{'ORF_'.$crt_NT_entries{$entry}{'protein_id'}}{$crt_NT_entries{$entry}{'protein_id'}.'_'.$crt_NT_entries{$entry}{'product'}}{$start_end_orientation_by_PID{$crt_NT_entries{$entry}{'protein_id'}}} = $ORF_names{$start_end_orientation_by_PID{$crt_NT_entries{$entry}{'protein_id'}}}; #generate the plot entry

												
											} else {
												#~ print "\t\t\t  $entry not in relation to previous results\n";
											}
										}	
									}
								} else {
									undef $crt_NT_ACC;
								}
							} 
							unless ($crt_NT_ACC) {
								$crt_NT_ACC = "origin_of_$crt_PID";
								$sequences{$sample_ID.'|||'.$sequence_header}{$crt_NT_ACC}{'NT'}{$crt_NT_ACC.'_'.$crt_PID_entries{'main'}{'definition'}}{'1_100'}= "could not reliably determine a nucleotide accession number for $crt_PID in $crt_PID_entries{'main'}{'source-db'}\nTaxonomy: ".$crt_PID_entries{'main'}{'taxonomy'};

								push (@associated_PIDs , $crt_PID);
								#~ print "\t\t\tcould not reliably determine a nucleotide accession number for $crt_PID in $crt_PID_entries{'main'}{'source-db'} . Will use $crt_NT_ACC as substitute!\n";
								my $crt_reference_ORF = $crt_NT_ACC.'_ORF_'.$crt_PID;
								#~ print "setting refernence for $sequence_header.'_rainbow_'.$crt_reference_ORF\n";
								$ORF_needs_a_rainbow{$sequence_header.'_rainbow_'.$crt_reference_ORF} = $crt_ORF;
								$reference_ORFs_AA{$crt_reference_ORF} = $crt_PID_entries{'main'}{'sequence'};
								$start_end_orientation_by_PID{$crt_PID}= '1_'. length($crt_PID_entries{'main'}{'sequence'}) *3 .'_forward';
								print "\t\t\t  running rainbow on $crt_reference_ORF vs $crt_ORF\n" ;
								#~ print "\t\t\t\tthat would be:\n\t\t\t\t$crt_ORF: $suspicious_ORFs_AA{$crt_ORF}\n";
								#~ print "\t\t\t\tvs.\n\t\t\t\t$crt_reference_ORF: $reference_ORFs_AA{$crt_reference_ORF}\n";
								$sequences{$sample_ID.'|||'.$sequence_header}{$crt_NT_ACC}{'ORF_'.$crt_PID}{$crt_PID.'_'.$crt_PID_entries{'main'}{'definition'}}{$start_end_orientation_by_PID{$crt_PID}} = $crt_PID_entries{'main'}{'definition'}; #generate the plot entry
							}
						}
						if ($crt_PID_entries{'main'}{'product'}){
							push(@{$potential_annotations{$crt_ORF}}, $crt_PID_entries{'main'}{'product'}) unless ($crt_PID_entries{'main'}{'product'} ~~ @{$potential_annotations{$crt_ORF}});
						}
					}
					#~ print "\n";
				} else {
					#~ print "\t\tno significant match yielded by any method\n\n";
				}
				my $joined_potential_annotations = "potential annotations:\n";
				if ($potential_annotations{$crt_ORF}){
					$joined_potential_annotations .= join("\n", sort (@{$potential_annotations{$crt_ORF}}));
					$joined_potential_annotations .=  "\n\nfound by: \n$crt_ORF_found_by_methods\n";
					#~ print "#potential annotations:\n$joined_potential_annotations\n#\n";
				} else {
					#~ print "#no potential annotations available\n#\n";
					$joined_potential_annotations = "no potential annotations available";
				}
				
				$crt_ORF =~ m/(.+)_(ORF_\d+)/;
				my $crt_ORF_number  = $2;
				
				$sequences{$sample_ID.'|||'.$sequence_header}{$sequence_header}{$crt_ORF_number}{$crt_ORF_number}{$crt_ORF_start_end_orientation} = "$joined_potential_annotations\n "; #generating plot entry for base NT sequence

				
			}
		}
		my $sequence_organization_details_file = File::Spec -> catfile ($TCC{'result_dir'}, $sample_ID.'_sequence_organization_details.csv'); 
		open (my $SEQUENCE_ORGANIZATION_DETAILS_FH, '>',  $sequence_organization_details_file) or die "could not write to '$sequence_organization_details_file': $!\n";
		###########creating the plots
		my $number_of_seqs; #will count the number of sequences in that plot
		my $max_length = 0; #the longest sequence will be stored to determine the needed width of the SVG
		my $rows_needed; #depending on the number of sequences and
		my $SVG_content; #will contain the SVG code
		my $row_start = 80; #start row for plot 
		my $rainbow_details;
		my $row_offset = $row_start; #the starting row for each sequence (y axis value)
		my $groupcount; #counts the number of groups (= number of suspicious sequences)
		my %color_blocks; #will contain the ranges for color blocks and the associated color colored by rainbow colors vong colors her 
					      #ranges will be calculated to be exact SVG x-values from the start of the sequence on
		my %color_blocks_reverse; #the same thing as above....but reverse
		my %color_block_order; #will store the order of ranges to loop through for each ORF
		my %ORF_plot_orientation;#will store the orientation in the plot for each ORF
		
		
		$SVG_content .= qq(<text x="100" y="40" font-weight="bold" font-size="2em" >$sample_ID<title>$crt_sample_info\n</title></text>\n); #adding the name of the group to the plot
		print "creating sequence organization plot for $sample_ID\n";
		foreach my $group (sort keys %sequences){
			(my $group_leader = $group) =~s/^[^\|]+\|{3}//;
			my @sequence_order = $group_leader;
			foreach my $sequence (sort keys %{$sequences{$group}}){
				push (@sequence_order, $sequence) unless ($sequence ~~ @sequence_order);
			}
			print "adding group $group to plot\n";
			$groupcount++;
			$SVG_content .= qq(<text x="100" y="$row_offset" font-weight="bold" font-size="1.5em" >$group_leader</text>\n); #adding the name of the group to the plot
			print $SEQUENCE_ORGANIZATION_DETAILS_FH "$sample_ID,$group_leader\n";
			$row_offset += 40; #next group...so we push the stuff further down so that is separated from the previous group
			my @group_leader_ORFs;
			foreach my $sequence (@sequence_order){
				print "\tadding sequence $sequence to plot\n"; 
				my %plot_rows; #key row number, value is an array containing $start_end_orientation
				my $rowcount = 2; #we start with 2 rows: 1st row is for the NT sequence (blue) and the 2nd row for the first ORF
				foreach my $sequence_type (sort keys %{$sequences{$group}{$sequence}}){
					foreach my $sequence_description (sort keys %{$sequences{$group}{$sequence}{$sequence_type}}){
						my @start_end_orientation = keys %{$sequences{$group}{$sequence}{$sequence_type}{$sequence_description}}; #that is: from_to_forward/reverse
						my @limits = split ('_', $start_end_orientation[0]);
						my $xstart = 100 + $limits[0]; #+100 to leave whitespace to the left
						my $width = $limits[1] -  $limits[0];
						my $yoffset;
						if ($sequence_type eq 'NT'){ #if the sequence is NT, the sequence will be represented as a blue bar
							$yoffset= $row_offset + 10;
							$SVG_content .= qq(<rect x="$xstart" y="$yoffset"  width="$width" height="15" style="fill:darkblue"><title>$sequences{$group}{$sequence}{$sequence_type}{$sequence_description}{$start_end_orientation[0]}\n</title></rect>\n);
							
							my $textyoffset = $yoffset - 5;
							$SVG_content .= qq(<text x="$xstart" y="$textyoffset" >$sequence_description</text>\n);
							$max_length = $limits[1] if ($limits[1] > $max_length); #to find out the length of the longest sequence for setting witdh limit of the plot
						} else {#this applies to all not NT.....i. e. ORFS
							#~ print "\t\tadding $sequence_description to plot\n";
							#most important step: try to puzzle the ORFs into the plot without them overlapping and also not wasting space by putting each ORF in a single row
							my $does_not_fit ; #will contain the information whether the ORF fits in the current row
							my $frame = $limits[2];
							for (my $plot_row = 2; $plot_row <= $rowcount;  $plot_row++){
								undef $does_not_fit;
								$yoffset = $row_offset + $plot_row * 35;
								unless ( $plot_rows{$plot_row}){ #if the row is empty, just put the ORF-range in the row
									push (@{$plot_rows{$plot_row}}, $start_end_orientation[0]);
									$does_not_fit = 'false';
								} else {
									foreach my $occupied (@{$plot_rows{$plot_row}}){
										my @occupied_limits = split ('_', $occupied);
										if (($limits[0] >= $occupied_limits[0]) and ($limits[1] <= $occupied_limits[1])){
											$does_not_fit = 'true';
											last;
										} elsif (($limits[0] <= $occupied_limits[0]) and ($limits[1] >= $occupied_limits[1])){
											$does_not_fit = 'true';
											last;
										} elsif (($limits[0]  + 5 < $occupied_limits[0]) and ($limits[1] + 5 < $occupied_limits[0])){
											$does_not_fit = 'false';
										} elsif (($limits[0] - 5 > $occupied_limits[1]) and ($limits[1] - 5 > $occupied_limits[1])){
											$does_not_fit = 'false';
										}  else {
											$does_not_fit = 'true';
											last;
										}
										
									}
									if ($does_not_fit eq 'false'){
										push (@{$plot_rows{$plot_row}}, $start_end_orientation[0]);
										last;
									}
								}
							}
							if ($does_not_fit eq 'true'){
								$rowcount++;
								$yoffset += 35;
								push (@{$plot_rows{$rowcount}}, $start_end_orientation[0]);
							}
							my $color;
							if ($frame eq 'forward'){
								$color = 'darkgrey';
							} else {
								$color = 'lightgrey';
							}
							
							my $textyoffset = $yoffset - 5;
							$SVG_content .= qq(<text x="$xstart" y="$textyoffset" >$sequence_description</text>\n);
							
							
							if ( $sequence =~ m/$group_leader/){
								print "\t\tpushing $sequence\_$sequence_description to group leader ORFs\n";
								if ($ORF_needs_a_rainbow{$sequence.'_'.$sequence_description}){
									$SVG_content .= qq(<rect x="$xstart" y="$yoffset" width="$width" height="10" style="fill:black"></rect>\n);
									#~ print "\t\tit will become the base for a rainbow!\n";
									#~ print "\t\tsequence is $suspicious_ORFs_AA{$sequence.'_'.$sequence_description}\n";
									my $crt_length = length($suspicious_ORFs_AA{$sequence.'_'.$sequence_description}) * 3;
									#~ print "\t\tsequence length in NT is $crt_length\n";
									my @color_block_ranges = 'NIL';
									my @colors = 'NIL';
									for (my $perc = 0; $perc <=100 - $color_resolution; $perc += $color_resolution){
										my $color_index = int($perc/100*244)+401;
										my $perc_plus_increment = $perc + $color_resolution; 
										my $from = int($perc/100*$crt_length) +1;
										my $to = int($perc_plus_increment/100*$crt_length);
										push(@color_block_ranges, $from.'_'.$to);
										push(@colors, $rainbow{$color_index});
									}	
									for (my $i = 1; $i < scalar @color_block_ranges; $i++){
										$color_blocks{$sequence.'_'.$sequence_description}{$color_block_ranges[$i]} = $colors[$i];
										$color_blocks_reverse{$sequence.'_'.$sequence_description}{$color_block_ranges[$i]} = $colors[-$i];
										#~ print "\t\t\t$color_block_ranges[$i] forward: $colors[$i], reverse: $colors[-$i]\n";
										push (@{$color_block_order{$sequence.'_'.$sequence_description}}, $color_block_ranges[$i]);
									}
									
									my %used_color_blocks;
									if ($frame eq 'forward'){
										%used_color_blocks = %{$color_blocks{$sequence.'_'.$sequence_description}};
										$ORF_plot_orientation{$sequence.'_'.$sequence_description} = 'forward';
									} else {
										%used_color_blocks = %{$color_blocks_reverse{$sequence.'_'.$sequence_description}};
										$ORF_plot_orientation{$sequence.'_'.$sequence_description} = 'reverse';
									}
									
									my $block_counter;
									foreach my $crt_color_block (@{$color_block_order{$sequence.'_'.$sequence_description}}){
										$block_counter++;
										#~ print "\t\t\t\t$crt_color_block:\t\t$used_color_blocks{$crt_color_block}\n";
										my @color_block_range = split ('_', $crt_color_block);
										my $color_block_xstart = $xstart + $color_block_range[0] -1;
										my $color_block_width = $color_block_range[1] - $color_block_range[0] + 1;
										$color_block_width-- if ($block_counter == scalar @{$color_block_order{$sequence.'_'.$sequence_description}});
										my $color_block_yoffset = $yoffset + 1;
										$SVG_content .= qq(<rect x="$color_block_xstart" y="$color_block_yoffset" width="$color_block_width" height="8" style="fill:white"></rect>\n);
										$SVG_content .= qq(<rect x="$color_block_xstart" y="$color_block_yoffset" width="$color_block_width" height="8" style="fill:$used_color_blocks{$crt_color_block}"></rect>\n);
									} 
									$SVG_content .= qq(<rect x="$xstart" y="$yoffset" width="$width" height="10" style="fill:white; opacity:0"><title>$sequences{$group}{$sequence}{$sequence_type}{$sequence_description}{$start_end_orientation[0]}</title></rect>\n);
									my @mod_seq_annotations = split ('found by:',$sequences{$group}{$sequence}{$sequence_type}{$sequence_description}{$start_end_orientation[0]});
									my $mod_seq_annotation;
									if ($mod_seq_annotations[1]){
										$mod_seq_annotations[0] =~ s/\n+/&/g;
										$mod_seq_annotations[0] =~ s/^\s+|\s+$//g;
										$mod_seq_annotations[0] =~ s/^&|&$//g;
										$mod_seq_annotations[1] =~ s/\n+/&/g;
										$mod_seq_annotations[1] =~ s/^\s+|\s+$//g;
										$mod_seq_annotations[1] =~ s/^&|&$//g;
										$mod_seq_annotations[0] =~ s/^potential annotations\:\&//;
										$mod_seq_annotation = "$mod_seq_annotations[0],$mod_seq_annotations[1]";
									} else {
										$mod_seq_annotations[0] =~ s/\n+/&/g;
										$mod_seq_annotations[0] =~ s/^\s|\s$//g;
										$mod_seq_annotations[0] =~ s/^&|&$//g;
										$mod_seq_annotation = "$mod_seq_annotations[0],NA";
									}
									print $SEQUENCE_ORGANIZATION_DETAILS_FH "$sample_ID,$sequence,$sequence_description,$mod_seq_annotation\n";

								}else {
										$SVG_content .= qq(<rect x="$xstart" y="$yoffset" width="$width" height="10" style="fill:$color"><title>$sequences{$group}{$sequence}{$sequence_type}{$sequence_description}{$start_end_orientation[0]}</title></rect>\n);
									my @mod_seq_annotations = split ('found by:',$sequences{$group}{$sequence}{$sequence_type}{$sequence_description}{$start_end_orientation[0]});
									my $mod_seq_annotation;
									if ($mod_seq_annotations[1]){
										$mod_seq_annotations[0] =~ s/\n+/&/g;
										$mod_seq_annotations[0] =~ s/^\s+|\s+$//g;
										$mod_seq_annotations[0] =~ s/^&|&$//g;
										$mod_seq_annotations[1] =~ s/\n+/&/g;
										$mod_seq_annotations[1] =~ s/^\s+|\s+$//g;
										$mod_seq_annotations[1] =~ s/^&|&$//g;
										$mod_seq_annotations[0] =~ s/^potential annotations\:\&//;
										$mod_seq_annotation = "$mod_seq_annotations[0],$mod_seq_annotations[1]";
									} else {
										$mod_seq_annotations[0] =~ s/\n+/&/g;
										$mod_seq_annotations[0] =~ s/^\s|\s$//g;
										$mod_seq_annotations[0] =~ s/^&|&$//g;
										$mod_seq_annotation = "$mod_seq_annotations[0],NA";
									}
									print $SEQUENCE_ORGANIZATION_DETAILS_FH "$sample_ID,$sequence,$sequence_description,$mod_seq_annotation\n";
									}
							} else {
								print "#####checking ", $sequence.'_'.$sequence_type, " for rainbow\n";
								#~ print "reference for: $group_leader.'_rainbow_'.$sequence.'_'.$sequence_type\n";
								if ($ORF_needs_a_rainbow{$group_leader.'_rainbow_'.$sequence.'_'.$sequence_type}){
									
									my $associated_ORF = $ORF_needs_a_rainbow{$group_leader.'_rainbow_'.$sequence.'_'.$sequence_type};
									if ($color_blocks{$associated_ORF}){
										my $associated_ORF_sequence = $suspicious_ORFs_AA{$associated_ORF};
										my $associated_ORF_sequence_reverse = reverse $associated_ORF_sequence;
										#~ print "\t\tit will get a rainbow based on $associated_ORF!\n";
										#~ sleep(1); 
										open (my $RAINBOW_DB, '>', 'rainbow_db_tmp.fasta') or die "could not write to 'rainbow_db_tmp.fasta': $!\n";
										print $RAINBOW_DB ">$associated_ORF\n$associated_ORF_sequence\n>$associated_ORF\_reverse\n$associated_ORF_sequence_reverse\n";
										close ($RAINBOW_DB);
										#~ $sequences{$crt_header.'_reversed'} = reverse $sequences{$crt_header};
										open (my $RAINBOW_QRY, '>', 'rainbow_qry_tmp.fasta') or die "could not write to 'rainbow_qry_tmp.fasta': $!\n";
										print $RAINBOW_QRY ">$sequence\_$sequence_type\n", $reference_ORFs_AA{$sequence.'_'.$sequence_type},"\n";
										close ($RAINBOW_QRY);
										#~ print "\t\t\tblasting...\n";
										#~ system(qq($TCC{'blastp'} -query $crt_sample_suspicious_fasta -out $crt_sample_suspicious_blasted $TCC{'blastp_settings'} -db $TCC{'blastp_db'} -num_threads $TCC{'nCPU'} -outfmt \"10 qseqid qlen qframe sframe length stitle pident nident qcovhsp qstart qend sstart send qseq sseq evalue score bitscore sacc\")) and die "Fatal: blast failed for '$crt_sample_suspicious_fasta': $!\n";
										system(qq($TCC{'makeblastdb'} -in rainbow_db_tmp.fasta -dbtype prot -out rainbow_db -title rainbow_db )) and die "Fatal: could not build blast database from 'rainbow_db'\: $!\n";
										#~ system(qq($TCC{'blastp'} -query rainbow_qry.tmp -out rainbow_results.csv -db rainbow_db -evalue $color_sensitivity -num_threads $TCC{'nCPU'} -outfmt \"10 qseqid qlen slen length stitle pident nident qcovhsp qstart qend sstart send qseq sseq evalue score bitscore sacc\")) and die "Fatal: blast failed for 'rainbow_qry.tmp ': $!\n";
										system(qq($TCC{'blastp'} -query rainbow_qry_tmp.fasta -out rainbow_results.csv -db rainbow_db -num_threads $TCC{'nCPU'} -outfmt \"10 qseqid qlen slen length stitle pident nident qcovhsp qstart qend sstart send qseq sseq evalue score bitscore sacc\")) and die "Fatal: blast failed for 'rainbow_qry.tmp ': $!\n";
										
										if (-s 'rainbow_results_sorted.csv') {
											system(qq(rm rainbow_results_sorted.csv));
										}

										#~ print "\t\t\tblast done\n";
										&sort_CSV(\'rainbow_results.csv', \'3', \'true');
										if (-s 'rainbow_results_sorted.csv'){
											$SVG_content .= qq(<rect x="$xstart" y="$yoffset" width="$width" height="10" style="fill:$color"></rect>\n);
											#~ print "\t\t\treading BLAST output...\n";
											open (my $BLAST_RESULTS, '<', 'rainbow_results_sorted.csv') or die "could not read from 'rainbow_results_sorted.csv': $!\n";
											my $crt_seq = 'Not initialized';
											my $identities;
											my $blast_lines;
											my $blast_areas;
											my $detailed_blastresults;
											while (my $line = <$BLAST_RESULTS>){
												chomp $line;
												print "$line,\n";
												$rainbow_details .= $line."\n";
												my @cols = split(',', $line);
												$cols[4] =~ m/(ORF_\d+)/;
												#~ die "here we have $cols[5]\n" if (int($cols[5]) eq 100);
												#~ print "\t\t\tdetected $cols[4] at $cols[5]% identity from $cols[10] to $cols[11] on current ORF from $cols[8] to $cols[9]\n";
												$detailed_blastresults .= "$cols[5]% identity from $cols[10] to $cols[11] vs. $1 from $cols[8] to $cols[9]\n";
												$identities += $cols[5];
												$blast_lines++;
												my $sstart = $cols[10]*3-2;
												my $qstart= $cols[8]*3-2;
												my $send= $cols[11]*3;
												my $qend = $cols[9]*3;
												my $qlen = $cols[1]*3;
												my $slen = $cols[2]*3;
												#~ my $plot_direction;
												my $hit_plot_start;
												if ($frame eq 'forward') {
													$hit_plot_start = $xstart +$qstart;
														
												} else {
													$hit_plot_start = $xstart + $qlen - $qend;
												}
												my $hit_plot_width =$qend - $qstart;
												my %used_color_blocks;
												if ($ORF_plot_orientation{$associated_ORF} eq 'forward'){
													if ($frame eq 'forward') {
														if ($cols[4] !~ m/reverse/) {
															%used_color_blocks = %{$color_blocks{$associated_ORF}};
														} else {
															%used_color_blocks = %{$color_blocks_reverse{$associated_ORF}};
														}
													 } else {
														 if ($cols[4] !~ m/reverse/) {
															%used_color_blocks = %{$color_blocks_reverse{$associated_ORF}};
														} else {
															%used_color_blocks = %{$color_blocks{$associated_ORF}};
														}
													 }
												 } else {
													if ($frame eq 'forward') {	
														if ($cols[4] !~ m/reverse/) {
															%used_color_blocks = %{$color_blocks{$associated_ORF}};
														} else {
															%used_color_blocks = %{$color_blocks_reverse{$associated_ORF}};
														}
													} else {	 
														if ($cols[4] !~ m/reverse/) {
															%used_color_blocks = %{$color_blocks{$associated_ORF}};
														} else {
															%used_color_blocks = %{$color_blocks_reverse{$associated_ORF}};
														}
													}
												}
												my %color_at_position;
												###what colors are in the matched part?
												for (my $i = $sstart; $i <= $send; $i++){
													#~ #qprint "checking between $sstart and $send\n";
													foreach my $color_block (@{$color_block_order{$associated_ORF}}){# color blocks of the original sequence
														my @start_end = split('_', $color_block);
														#~ #print "\tchecking between $start_end[0] and $start_end[1]\n";
														if (($i >= $start_end[0]) and ($i <= $start_end[1])){
																$color_at_position{$i} = $used_color_blocks{$color_block};
															last;
														}
													}
												}
												
												
												my $crt_color;
												#~ $blast_areas .=  "$sstart\_to_$send\_is_$qstart\_to_$qend\____";
												my $opacity=sprintf ("%.2f",$cols[5]/ 100);
												#~ my $pinkyoffset = $yoffset -5;
												#~ $SVG_content .= qq(<rect x="$hit_plot_start" y="$pinkyoffset" width="$hit_plot_width" height="20" style="fill:white; opacity:0.1; stroke-width:1; stroke:black"><title>blast_area of result $blast_lines: $sstart to $send is $qstart to $qend in plot direction $frame\nbased on  $cols[10] to $cols[11] is $cols[8] to $cols[9] vs $cols[4]:\n\n$detailed_blastresults </title></rect>\n);
												
												my $color_block_xstart;
												my $color_block_yoffset= $yoffset + 1;
												my $position_counter;
												foreach my $position (sort {$a <=> $b} keys %color_at_position){
													if ($crt_color){
														if ($crt_color ne $color_at_position{$position}){
															$SVG_content .= qq(<rect x="$color_block_xstart" y="$color_block_yoffset" width="$position_counter" height="8" style="fill:white"></rect>\n);
															$SVG_content .= qq(<rect x="$color_block_xstart" y="$color_block_yoffset" width="$position_counter" height="8" style="fill:$crt_color; opacity:$opacity"></rect>\n);
															$color_block_xstart = $position + $hit_plot_start - $sstart;
															undef ($position_counter);
															$crt_color = $color_at_position{$position};
														}
													} else {
														$color_block_xstart = $position + $hit_plot_start - $sstart;
														$crt_color = $color_at_position{$position};
													}
													$position_counter++;
												}
												
												$SVG_content .= qq(<rect x="$color_block_xstart" y="$color_block_yoffset" width="$position_counter" height="8" style="fill:white"></rect>\n);
												$SVG_content .= qq(<rect x="$color_block_xstart" y="$color_block_yoffset" width="$position_counter" height="8" style="fill:$crt_color; opacity:$opacity"></rect>\n);
												$SVG_content .= qq(<rect x="$hit_plot_start" y="$yoffset" width="$hit_plot_width" height="10" style="fill:pink; opacity:0"><title>$cols[5]% identity from $cols[10] to $cols[11] vs. $1 from $cols[8] to $cols[9]\n</title></rect>\n);

	
											}
											close ($BLAST_RESULTS);
											
											#~ my $identity_average = int($identities/$blast_lines);
											#~ my $blastyoffset = $textyoffset +30;
											#~ $SVG_content .= qq(<rect x="$xstart" y="$yoffset" width="$width" height="10" style="fill:white; opacity:0"><title>$identity_average% average identity to $associated_ORF:\n$sequences{$group}{$sequence}{$sequence_type}{$sequence_description}{$start_end_orientation[0]}:\n\n$detailed_blastresults </title></rect>\n);
											
										} else {
											print "#####could not find blast results for $associated_ORF!\n";
											$SVG_content .= qq(<rect x="$xstart" y="$yoffset" width="$width" height="10" style="fill:$color"><title>$sequences{$group}{$sequence}{$sequence_type}{$sequence_description}{$start_end_orientation[0]}</title></rect>\n);
										}
										
									} else {
										print "#####could not find stuff for $associated_ORF!\n";
										$SVG_content .= qq(<rect x="$xstart" y="$yoffset" width="$width" height="10" style="fill:$color"><title>$sequences{$group}{$sequence}{$sequence_type}{$sequence_description}{$start_end_orientation[0]}</title></rect>\n);
									}
								} else {
									$SVG_content .= qq(<rect x="$xstart" y="$yoffset" width="$width" height="10" style="fill:$color"><title>$sequences{$group}{$sequence}{$sequence_type}{$sequence_description}{$start_end_orientation[0]}</title></rect>\n);
								}
							}
						}
					}
				}
				$rows_needed += $rowcount + 1;
				$number_of_seqs++;
				$row_offset += $rowcount * 35+ 80;
			}
			$row_offset += 50;
			
			undef %color_blocks;  
			undef  %color_blocks_reverse; 
			undef  %color_block_order; 
			undef %ORF_plot_orientation;
			
			
		}
		print "will create an SVG containing $number_of_seqs sequences in $rows_needed rows\n";
		my $width = $max_length + 300;
		
		my $height = $rows_needed * 35 + $row_start + $number_of_seqs * 70 + $groupcount * 40;
		my $SVG_header = qq(<svg xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.w3.org/2000/svg" height="$height" width="$width" version="1.1" xmlns:cc="http://creativecommons.org/ns#" xmlns:dc="http://purl.org/dc/elements/1.1/">);
		my $SVG_end = qq(</svg>);
		my $SVG_file = File::Spec -> catfile ($TCC{'result_dir'}, $sample_ID.'_sequence_organization.svg'); 
		open (my $SVG_FH, '>', $SVG_file) or die "could not write to '$SVG_file': $!\n";
		print $SVG_FH $SVG_header, "\n";
		print $SVG_FH $SVG_content;
		print $SVG_FH $SVG_end;
		print $SCAVENGED_SAMPLES $sample_ID, "\n";
		close ($SVG_FH);
		
		
		
		my $rainbow_details_file = File::Spec -> catfile ($TCC{'result_dir'}, $sample_ID.'_rainbow_details.csv'); 
		open (my $RAINBOW_DETAILS_FH, '>', $rainbow_details_file) or die "could not write to '$rainbow_details_file': $!\n";
		print $RAINBOW_DETAILS_FH $rainbow_details, "\n";
		close ($RAINBOW_DETAILS_FH);
		close ($SEQUENCE_ORGANIZATION_DETAILS_FH);
		
		my $t2 = time();
		my $tdiff = $t2 - $t1;
		print "needed $tdiff seconds to scavenge\n";
		
	}
} else {
	print 'all ', scalar keys %infection_status, " samples seem to be clean!\n";
}


(my $evaluated_sample_library = $sample_library) =~s/\.csv/_evaluated.csv/;
open (my $EVALUATED_SAMPLE_LIBRARY, '>', $evaluated_sample_library ) or die "could not write to'$evaluated_sample_library ': $!\n";
my $joined_header = join(',', @{$sample_library_entries{'header'}});
print $EVALUATED_SAMPLE_LIBRARY "$joined_header\n";

foreach my $sample (keys %infection_status){
	if ($infection_status{$sample} eq 'clean'){
		push (@{$sample_library_entries{$sample}}, 0);
	}
	
	my $joined_info = join(',', @{$sample_library_entries{$sample}});
	print $EVALUATED_SAMPLE_LIBRARY "$joined_info\n";
	
}


print "TRAVIS Scavenger completed\n";


###subroutines

sub fetch_gbx {
	#download NCBI efetch data in genebank XML format for a given NT accession number and parse it for associated PIDs
	#PIDs will be stored in an array
	my $sref_accession_number=$_[0]; #a NT accession number
	my $href_crt_entries = $_[1]; #in: undefined, out: defined
	my $sref_database = $_[2];
	my $sref_local_database = $_[3];
	print "\t\tfetching features for $$sref_accession_number from NCBI...\n";
	#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001542&rettype=gb&retmode=xml
	#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=NP_056795&rettype=gb&retmode=xml
	my $gb_xml = File::Spec -> catfile ($$sref_local_database, $$sref_accession_number.'_'.$$sref_database.'_genebank.xml');
	my $gb_entry;
	if (-s $gb_xml){
		#~ print "\t\treading existing genebank xml...\n";
		open (my $GB_XML, '<', $gb_xml) or die "could not read from '$gb_xml': $!\n";
		while (my $line = <$GB_XML>){
			$gb_entry .= $line;
		}
		close ($GB_XML);		
	} else {
		print "\t\tretrieving genebank xml from NCBI...\n";
		
		my $test;
		my $curl_tmp = 'curl_thread3.tmp';
		for (my $try = 1; $try <= 10; $try++){
			print "\t\t\tcurl try $try...\n";
			system(qq(curl -s --retry 10 --retry-delay 5  -o $curl_tmp -k -L "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$$sref_database&id=$$sref_accession_number&rettype=gb&retmode=xml" )) and print "curl failed at '$$sref_accession_number': $!\n";
			$test = qx/grep 'GBQualifier_name' $curl_tmp/;
			#~ print "test: $test\n";
			sleep(1);
			last if ($test =~ m/GBQualifier_name/);
		
		}
		
		undef $test if ($test !~ m/GBQualifier_name/);
		
		sleep(1);
		if ((-s $curl_tmp) and ($test)){
			#~ print "reading $curl_tmp\n";
			&slurp_file(\$curl_tmp,\$gb_entry);
			system(qq(mv $curl_tmp $gb_xml));
		} else {
			system (qq(rm $curl_tmp)) if (-s $curl_tmp);
			print "!!!!!!!!!!!!!could not fetch features for $$sref_accession_number from NCBI...\n";
			sleep(1); 
		}
	}
	
	
	
	
	if (!$gb_entry){
		sleep(1); 
		$href_crt_entries = 'failed';
	} else {
		my $entry_count;
		my %entry_elements;
		my @gb_entry_lines= split("\n",$gb_entry); #split the entry to read it line by line
		for (my $i=1; $i < scalar @gb_entry_lines; $i++){
			next if ($gb_entry_lines[$i] !~ m/GB/);
			if ($gb_entry_lines[$i] =~ m/<GBFeature>/){
			} elsif ($gb_entry_lines[$i] =~ m/<\/GBFeature>/){
					$entry_count++;
					%{$href_crt_entries -> {$entry_count}}= %entry_elements;
				undef %entry_elements;
			} elsif ($gb_entry_lines[$i] =~ m/<GBFeature_key>(.*)<\/GBFeature_key>/){
				#~ print "\tfound feature key: $1\n";
				$entry_elements{'key'} = $1;
			} elsif ($gb_entry_lines[$i] =~ m/<GBInterval_from>[^\d]*(\d+)[^\d]*<\/GBInterval_from>/){
				#~ print "\tfound feature from: $1\n";
				$entry_elements{'from'} = $1;
			} elsif ($gb_entry_lines[$i] =~ m/<GBInterval_to>[^\d]*(\d+)[^\d]*<\/GBInterval_to>/){
				#~ print "\tfound feature from: $1\n";
				$entry_elements{'to'} = $1;
			} elsif ($gb_entry_lines[$i] =~ m/<GBQualifier_name>(.*)<\/GBQualifier_name>/){
				my $name = $1;
				$i++;
				$gb_entry_lines[$i] =~ m/<GBQualifier_value>(.*)<\/GBQualifier_value>/;
				if ($1) {
					$entry_elements{$name} = $1;
					if ($name =~ m/host/i) {
						#~ print "$entry_elements{$name} !\n";
						$href_crt_entries -> {'main'}{'host'}= $entry_elements{$name} ;
					}
				} else {
					$entry_elements{$name} = 'NA';
				}
				
					
			#<GBSeq_source-db>REFSEQ: accession NC_001542.1</GBSeq_source-db>
			} elsif ($gb_entry_lines[$i] =~ m/<GBSeq_source-db>(.*)<\/GBSeq_source-db>/){
				$href_crt_entries -> {'main'}{'source-db'}= $1;
			} elsif ($gb_entry_lines[$i] =~ m/<GBQualifier_name>product<\/GBQualifier_name>/){
				$gb_entry_lines[$i+1] =~ m/<GBQualifier_value>(.*)<\/GBQualifier_value>/;
				$href_crt_entries -> {'main'}{'product'}= $1;
			} elsif ($gb_entry_lines[$i] =~ m/ <GBSeq_taxonomy>(.*)<\/GBSeq_taxonomy>/){
				$href_crt_entries -> {'main'}{'taxonomy'}= $1;
				#<GBSeq_primary-accession>AF013368</GBSeq_primary-accession>
			}elsif ($gb_entry_lines[$i] =~ m/ <GBSeq_primary-accession>(.*)<\/GBSeq_primary-accession>/){
				$href_crt_entries -> {'main'}{'primary-accession'}= $1;
			} elsif ($gb_entry_lines[$i] =~ m/ <GBSeq_sequence>(.*)<\/GBSeq_sequence>/){
				$href_crt_entries -> {'main'}{'sequence'}= uc ($1);
			} elsif ($gb_entry_lines[$i] =~ m/ <GBSeq_length>(.*)<\/GBSeq_length>/){
				$href_crt_entries -> {'main'}{'length'}= $1;
			} elsif ($gb_entry_lines[$i] =~ m/ <GBSeq_definition>(.*)<\/GBSeq_definition>/){
				($href_crt_entries -> {'main'}{'definition'}= $1) =~ s/[^a-zA-Z0-9_\-]/_/g;
				$href_crt_entries -> {'main'}{'definition'} =~ s/_$//;
			}
		}
		sleep(1); 
	}
}


sub  slurp_file{
	my $sref_file = $_[0];
	my $sref_container = $_[1];
	open (my $FILE, '<', $$sref_file) or print "could not read from '$$sref_file': $!\n";
	while (my $line = <$FILE>){
		$$sref_container  .= $line;
	} 
}



sub fetch_sequence{
	#download NCBI efetch data in fasta format for a given NT accession number and parse it for description and sequence
	my $sref_accession_number = $_[0]; #accession number or PID
	my $sref_crt_description = $_[1]; #in: defined or undefined, out: defined
	my $sref_crt_sequence = $_[2]; #in: defined or undefined, out: defined
	my $sref_database = $_[3]; #'nuccore' or 'protein'
#	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001542&rettype=fasta";

	my $test;
	my $curl_tmp = 'curl.tmp';
	for (my $try = 1; $try <= 10; $try++){
		#~ print "try $try...\n";
		system(qq(curl  --retry 10 --retry-delay 5 -s -o $curl_tmp -k -L "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$$sref_database&id=$$sref_accession_number&rettype=fasta")) and print "curl failed at '$$sref_accession_number': $!\n";
		$test = qx/grep '>' $curl_tmp/;
		#~ print "test: $test\n";
		sleep(1);
		last if ($test =~ m/>/);
		
	}
		
	undef $test if ($test !~ m/>/);
		
	sleep(1);
	if ((-s $curl_tmp) and ($test)){

		# print "reading $curl_tmp\n";
		undef $$sref_crt_sequence;
		my $gb_entry;
		&slurp_file(\$curl_tmp,\$gb_entry);
		my @lines = split("\n", $gb_entry);
		# print Dumper @lines;
		for (my $l=1; $l < scalar @lines;$l++){
			chomp $lines[$l];
			# print "$lines[$l]\n";
			$$sref_crt_sequence .= $lines[$l];
		}
		($$sref_crt_description = $lines[0]) =~ s/[^a-zA-Z0-9 _\-\[\]]|$$sref_accession_number[^\s]* |\s*\[[^\]]+\]//g;
		$$sref_crt_description =~ s/\s/_/g;
		system(qq(rm $curl_tmp));
	} else {
		print "failed to retrieve sequence for $$sref_accession_number\n";
		
	}
	sleep(1);  #NCBI does not like too many failed attempts in a short amount of time ;)
}

sub check_dir{
	#creates it if it does not exist
	my $sref_dir = $_[0]; #directory path
	unless (-d $$sref_dir){ #check if directory exists
		make_path($$sref_dir) or die "could not create '$$sref_dir': $!\n";
		print "'$$sref_dir' created\n";
	}	
}



sub read_TCC{
	#reads a TRAVIS Control Center file and stores settings in a hash
	my $sref_TCC_file = $_[0]; #file path
	my $href_TCC = $_[1]; #in: undefined, out: defined
	open (my $TCC_FH , '<', $$sref_TCC_file) or die "could not read from '$$sref_TCC_file': $!\n";
	while (my $line = <$TCC_FH>) {
		chomp $line;
		next if ($line =~ /^\#/); #skip comment lines
		next if ($line =~ /^\s*$/); #skip empty lines
		$line =~ s/"|'|//gi; #get rid of quotes 
		$line =~ s/\s+$//gi; #and trailing spaces 
		my @config_elements = split(',', $line); #split line into parameter and value
		#~ print "$config_elements[0]:$config_elements[1]\n";
		$href_TCC -> {$config_elements[0]}=$config_elements [1];
	}
	close $TCC_FH;
}

sub sort_CSV{
	my $sref_in_CSV = $_[0];
	my $sref_sort_by_column = $_[1];
	my $sref_ascending = $_[2];
	(my $out_CSV = $$sref_in_CSV)  =~ s/\.csv/_sorted.csv/;
	my %entries;
	my %values;
	#~ print "$$sref_sort_by_column\n";
	
	open (my $IN_CSV, '<', $$sref_in_CSV) or die "could not read from '$$sref_in_CSV': $1\n";
	open (my $OUT_CSV, '>', $out_CSV) or die "could not write to '$out_CSV': $1\n";
	#~ my $header;
	my $entrycount;
	while (my $line = <$IN_CSV>){
		chomp $line;
		#~ chomp $line =;
		#~ print $line, "\n";
		#~ if (!$header){
			#~ $header = $line;
			#~ print "$header \n";
			#~ print $OUT_CSV "$header \n";
			#~ next;
		#~ }
		$entrycount++;
		my @cols = split (',', $line);
		
		(my $key = $cols[$$sref_sort_by_column]) =~ s/\.\d+$//;
		$entries{'entry'.$entrycount} = $line;
		$values{'entry'.$entrycount} =  $key;
		
		#~ print "entry$entrycount at $cols[$$sref_sort_by_column] starting with $cols[0]\n";
	}
	#~ return;
	#~ print "read $entrycount entries\n";
	my @sorted_keys;
	if ($$sref_ascending eq 'true'){
		#~ sort { $planets{$a} <=> $planets{$b} } keys %planets)
		@sorted_keys = sort ( {$values{$a} <=> $values{$b}} keys %values);
	} else {
		@sorted_keys = sort ( {$values{$b} <=>  $values{$a}} keys %values);
	}
	
	my $lalala;
	foreach my $entry (@sorted_keys){
		$lalala++;
		
		#~ print "$lalala: $entry with $values{$entry}\n\t$entries{$entry}\n";
		print $OUT_CSV "$entries{$entry}\n";
	}
	close ($OUT_CSV);
	#~ print "wrote $lalala entries\n";
}

sub get_times {
	my $sref_done =$_[0];
	(my $out = $ $sref_done) =~ s/.csv/_times.csv/;
	my %times_by_method;
	my %times_by_sample;
	open (my $DONE, '<', $$sref_done) or die "could not read from '$$sref_done': $!\n";
	open (my $OUT, '>', $out) or die "could not write to'$out': $!\n";
	while (my $line = <$DONE>){
		chomp $line;
		my @cols = split (',', $line);
		next if ($times_by_sample{$cols[0].$cols[3].$cols[5]});
		next if ($cols[4] eq 'NA');
		$times_by_sample{$cols[0].$cols[3].$cols[5]} = $cols[4];
		$times_by_method{$cols[5]} += $cols[4];
	}
	print $OUT "method,seconds,minutes,hours,days\n";
	foreach my $method (keys %times_by_method){
		my $mins = $times_by_method{$method} / 60;
		my $hours = $mins / 60;
		my $days = $hours / 24;
		print  $OUT  "$method,$times_by_method{$method},$mins,$hours,$days\n";
	}
}
