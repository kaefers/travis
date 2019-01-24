#!/usr/bin/env perl
use warnings;
use strict;
use LWP::Simple;
use autodie;
use File::Path qw(make_path);	
use File::Basename;
use File::Spec;
use Data::Dumper;


unless (@ARGV){
	die "usage: perl TRAVIS_Core_v0.6a.pl <TRAVIS_control_center.csv>\n";
}
my $TCC_file = shift @ARGV;
#~ my $TCC_file = '/home/skaefer/reo_full/reo_full_TCC.csv';
#~ my $TCC_file = '/home/skaefer/sync/work/parkinson/synuclein_TCC.csv';

print '#' x 70,"\n";
print '#' x 70,"\n";
print "\t\t\tThis is TRAVIS Core v0.8a\n";
print '#' x 70,"\n";
print '#' x 70,"\n";

###reading TRAVIS Control Center
my %TCC; #contains configuration parameters from $TCC_file 
&read_TCC(
	\$TCC_file, #file path
	\%TCC #undef hash
);
my @supported_search_tools = ('blastp', 'hmmer', 'jackhmmer', 'mmseqs');
print "connected to TRAVIS Control Center\n";


my %status_of_ORF_fasta;
###loop through sample library and get the ORFs
print "checking ORF_dir..\n";
&check_dir(\$TCC{'ORF_dir'});
my @skipped; #collect filenames that are in the sample file but cannot be found
my $ORF_fastas_content; #collect the ORF fasta paths that can actually be searched
open (my $SAMPLE_LIBRARY, '<', $TCC{'sample_library'}) or die "could not read from $TCC{'sample_library'}: $!\n";

my $c1;
while (my $line = <$SAMPLE_LIBRARY>){
	chomp $line;
	next if ($line =~ m/^#/);
	
	my @cols = split (',', $line);
	my $sample_file = File::Spec -> catfile ($TCC{'sample_dir'} , $cols[0]);
	#~ my $crt_sample_ORF_data = File::Spec -> catfile ($TCC{'ORF_dir'} , $cols[1].'ORFs.fasta');
	unless (-s $sample_file){
		push (@skipped, $cols[0]);
		next;
	}
	
	#~ my $min_ORF_length = 50;

	#~ (my $crt_filename = basename($sample_file)) =~ s/\.[^\.]+$//;
	my $sample_ORF_data = File::Spec -> catfile($TCC{'ORF_dir'} , $cols[1].'_ORF_data.csv');
	unless ((-s $sample_ORF_data) or ($TCC{'purge_data'})){
		print "getting ORF data from $sample_ORF_data...\n";
		&get_sample_ORF_data(\$sample_file, \$sample_ORF_data, \'fasta', \$TCC{'min_ORF_length'},\$TCC{'max_ORF_length'});
		$c1++;
		#~ last if ($c1 >= 5 );
	} else {
		print "$sample_ORF_data already exists...skipping...\n";
	}
	my $sample_ORF_fasta = File::Spec -> catfile($TCC{'ORF_dir'} , $cols[1].'_ORFs_AA.fasta');
	my $sample_ORF_fasta_basename = basename($sample_ORF_fasta);
	unless (-s $sample_ORF_fasta ){
		if (-s $sample_ORF_data){
			print "converting $sample_ORF_data to $sample_ORF_fasta...\n";
			&ORF_data_2_fasta(\$sample_ORF_data, \$sample_ORF_fasta);
			$ORF_fastas_content .= $sample_ORF_fasta."\n";
		} else {
			print "no ORFs found that match the criteria\n";
		}
	} else {
		print "$sample_ORF_fasta already exists...skipping...\n";
		$ORF_fastas_content .= $sample_ORF_fasta_basename."\n";
		#~ next;
	}
	$status_of_ORF_fasta{$sample_ORF_fasta_basename} = 'clean';
	push (@{$status_of_ORF_fasta{'order'}}, $sample_ORF_fasta_basename);
	#~ print "read $sample_file\n";
	print '#' x 70,"\n";
}
###
print "checking result_dir..\n";
&check_dir(\$TCC{'result_dir'});	
my $ORF_fastas = File::Spec -> catfile ($TCC{'result_dir'}, 'ORF_fastas.txt');
&write_file(\$ORF_fastas ,\$ORF_fastas_content); 
foreach my $skipped (@skipped){
	print "skipped $skipped because it was empty or non existent!\n";
}



my $ORF_fastas_done = File::Spec -> catfile ($TCC{'result_dir'}, 'ORF_fastas_done.txt');
my $ORF_results = File::Spec -> catfile ($TCC{'result_dir'},  'ORF_results.csv');


my %already_done;
my $already_done_ORF_results;
if ((-s $ORF_fastas_done) and ($TCC{'resume_calculation'})){
	print "will try to resume calculation!\n";
	open (my $ORF_FASTAS_DONE, '<', $ORF_fastas_done) or die "could not read from '$ORF_fastas_done': $!\n";
	while (my $line = <$ORF_FASTAS_DONE>){
		chomp $line;
		print $line,   "\n";
		my @cols = split (',', $line);
		$already_done{$cols[0]} = $cols[1];
		$status_of_ORF_fasta{$cols[0]} = $cols[1];
	}
	close($ORF_FASTAS_DONE);
	open (my $ORF_RESULTS_DONE, '<',  $ORF_results) or die "could not read from '$ORF_results': $!\n";
	while (my $line = <$ORF_RESULTS_DONE>){
		chomp $line;
		#~ print $line,   "\n";
		my @cols = split (',', $line);
		my $crt_ORF_fasta = File::Spec -> catfile  ($TCC{'ORF_dir'}, $cols[0].'_ORFs_AA.fasta');
		$already_done_ORF_results .= $line."\n" if $already_done{$crt_ORF_fasta};
	}
	close($$ORF_RESULTS_DONE);
}

#~ print Dumper %already_done;
open (my $ORF_FASTAS_DONE, '>', $ORF_fastas_done) or die "could not write to '$ORF_fastas_done': $!\n";
open (my $ORF_RESULTS, '>', $ORF_results) or die "could not write to '$ORF_results': $!\n";
print $ORF_RESULTS $already_done_ORF_results if ( $already_done_ORF_results);

my $collection=  File::Spec -> catfile ('/home/skaefer/tmp/termites/references/testi', '*AA*collection*');
my @cfiles = glob $collection;
if (@cfiles){	
	print "removing previous collections...\n";
	system(qq(rm $collection));
} else {
	print "no previous collections found that should be deleted\n";
}

my %collection_contents;#dim1: all/main_positive, dim2: search tool, dim3: main/company, = array of fasta paths
my %reference_descriptions;


open (my $TTT_FH, '<', $TCC{'TTT'}) or die "could not read from '$TCC{'TTT'}': $!\n";
while (my $line = <$TTT_FH>){
	chomp $line;
	
	#~ print $line, "\n";
	my @cols = split (',', $line);
	next if ($cols[7] ne 'on');

	#~ print "preparing search for sequence(s) in $cols[0]\n";
	$reference_descriptions{$cols[0]} = $cols[1];
	$cols[1] =~ m/^([^_]+)/;
	my $main_or_company = $1;
	#~ my $main_or_company = 'main';
	#~ print "$main_or_company\n";
	my @search_tools = split ('&', $cols[6]);
	foreach my $search_tool (@search_tools){
		die "$search_tool is not supported!\n" unless ($search_tool ~~ @supported_search_tools);
		my $reference_file_path;
		if ($search_tool eq 'hmmer'){
			$reference_file_path = File::Spec -> catfile ($TCC{'reference_fastas'}, $cols[3]);
		}else {
			$reference_file_path = File::Spec -> catfile ($TCC{'reference_fastas'}, $cols[2]);
		}	
		push (@{$collection_contents{$cols[5]}{$search_tool}{$main_or_company}{$cols[0]}}, $reference_file_path);
	}
}




my $c;
foreach my $crt_sample (@{$status_of_ORF_fasta{'order'}}){
	print '#' x 70,"\n";
	print '#' x 70,"\n";
	print '#' x 70,"\n";
	print "searching in $crt_sample...\n";
	$c++;
	#~ last if ($c >= 7);
	if ($already_done{$crt_sample}){
		print "already did that! Skipping...\n";
		print $ORF_FASTAS_DONE "$crt_sample,$already_done{$crt_sample}\n";
		next;
	}
	my $crt_sample_filepath = File::Spec -> catfile ($TCC{'ORF_dir'}, $crt_sample);
	my @suspicious_sequences;
	foreach my $search_tool (sort keys %{$collection_contents{'all'}}){
	#~ foreach my $search_tool ('mmseqs'){
		print '#' x 70,"\n";
		print "\trunning $search_tool for main sequences...\n";
		if ($search_tool eq 'hmmer'){
			foreach my $reference (keys %{$collection_contents{'all'}{$search_tool}{'main'}}){
				foreach my $reference_file (@{$collection_contents{'all'}{$search_tool}{'main'}{$reference}}){
					print "\tvs. $reference_file\n";
					#~ print "\t\t\t###running hmmer...";
					&run_hmmer(
					\$reference, #name of the reference file
					\$reference_descriptions{$reference} , #description
					\$reference_file, #reference alignment
					$ORF_RESULTS, #filehandle
					\%TCC, #path to config
					\%status_of_ORF_fasta, #contains all ORF fastas as keys
					\$crt_sample_filepath,
					\@suspicious_sequences
					);
				}
			}
		} else {
			my $crt_reference_fasta = File::Spec -> catfile ($TCC{'reference_fastas'}, $search_tool.'_AA_all_main.fasta');
			print "\tjoining fastas to $crt_reference_fasta\n";
			unless (-s $crt_reference_fasta){ 
				&join_fastas(\%{$collection_contents{'all'}{$search_tool}{'main'}}, \$crt_reference_fasta);
			}
			if ($search_tool eq 'jackhmmer'){
				#~ print "\t\t\t###running jackhmmer...";
				&run_jackhmmer(
					\'jackhmmer_AA_all_main_collection', #name of the reference file
					\'jackhmmer_AA_all_main', #description
					\$crt_reference_fasta , #reference alignment
					$ORF_RESULTS, #filehandle
					\%TCC, #path to config
					\%status_of_ORF_fasta, #contains all ORF fastas as keys
					\$crt_sample_filepath,
					\@suspicious_sequences
				);	
			}elsif ($search_tool eq 'blastp'){
				#~ print "\t\t\t###running blastp...";
				&run_blastp(
					\'blastp_AA_all_main_collection', #name of the reference file
					\'blastp_AA_all_main', #description
					\$crt_reference_fasta , #reference alignment
					$ORF_RESULTS, #filehandle
					\%TCC, #path to config
					\%status_of_ORF_fasta, #contains all ORF fastas as keys
					\$crt_sample_filepath,
					\@suspicious_sequences
				);	
			}elsif ($search_tool eq 'mmseqs'){
				#~ print "\t\t\t###running mmseqs...";
				&run_mmseqs(
					\'mmseqs_AA_all_main_collection', #name of the reference file
					\'mmseqs_AA_all_main', #description
					\$crt_reference_fasta , #reference alignment
					$ORF_RESULTS, #filehandle
					\%TCC, #path to config
					\%status_of_ORF_fasta, #contains all ORF fastas as keys
					\$crt_sample_filepath,
					\@suspicious_sequences
				);	
			}
				
		}
		print '#' x 70,"\n";
		if ($search_tool eq 'hmmer'){
			foreach my $reference (keys %{$collection_contents{'all'}{$search_tool}{'company'}}){
				foreach my $reference_file (@{$collection_contents{'all'}{$search_tool}{'company'}{$reference}}){
					print "\tvs. $reference_file\n";
					print "\t\t\t###running hmmer...";
					&run_hmmer(
					\$reference, #name of the reference file
					\$reference_descriptions{$reference} , #description
					\$reference_file, #reference alignment
					$ORF_RESULTS, #filehandle
					\%TCC, #path to config
					\%status_of_ORF_fasta, #contains all ORF fastas as keys
					\$crt_sample_filepath,
					\@suspicious_sequences
					);
				}
			}
		} else {
			my $crt_reference_fasta = File::Spec -> catfile ($TCC{'reference_fastas'}, $search_tool.'_AA_all_company.fasta');
			print "\tjoining fastas to $crt_reference_fasta\n";
			unless (-s $crt_reference_fasta){ 
				&join_fastas(\%{$collection_contents{'all'}{$search_tool}{'company'}}, \$crt_reference_fasta);
			}
			if ($search_tool eq 'jackhmmer'){
				#~ print "\t\t\t###running jackhmmer...";
				&run_jackhmmer(
					\'jackhmmer_AA_all_company_collection', #name of the reference file
					\'jackhmmer_AA_all_company', #description
					\$crt_reference_fasta , #reference alignment
					$ORF_RESULTS, #filehandle
					\%TCC, #path to config
					\%status_of_ORF_fasta, #contains all ORF fastas as keys
					\$crt_sample_filepath,
					\@suspicious_sequences
				);	
			}elsif ($search_tool eq 'blastp'){
				#~ print "\t\t\t###running blastp...";
				&run_blastp(
					\'blastp_AA_all_company_collection', #name of the reference file
					\'blastp_AA_all_company', #description
					\$crt_reference_fasta , #reference alignment
					$ORF_RESULTS, #filehandle
					\%TCC, #path to config
					\%status_of_ORF_fasta, #contains all ORF fastas as keys
					\$crt_sample_filepath,
					\@suspicious_sequences
				);	
			}elsif ($search_tool eq 'mmseqs'){
				#~ print "\t\t\t###running mmseqs...";
				&run_mmseqs(
					\'mmseqs_AA_all_company_collection', #name of the reference file
					\'mmseqs_AA_all_company', #description
					\$crt_reference_fasta , #reference alignment
					$ORF_RESULTS, #filehandle
					\%TCC, #path to config
					\%status_of_ORF_fasta, #contains all ORF fastas as keys
					\$crt_sample_filepath,
					\@suspicious_sequences
				);	
			}
		}
		
		
	}
	if ($status_of_ORF_fasta{$crt_sample} eq 'clean'){
		print "no hits found for main set. Will skip company sequence search for main_positive\n";
		
	} else {
		foreach my $search_tool (keys %{$collection_contents{'main_positive'}}){
			print '#' x 70,"\n";
			print "\trunning $search_tool for company sequences in main_positive...\n";
			if ($search_tool eq 'hmmer'){
				foreach my $reference (keys %{$collection_contents{'main_positive'}{$search_tool}{'company'}}){
					foreach my $reference_file (@{$collection_contents{'main_positive'}{$search_tool}{'company'}{$reference}}){
						print "\tvs. $reference_file\n";
						#~ print "\t\t\t###running hmmer...";
						&run_hmmer(
						\$reference, #name of the reference file
						\$reference_descriptions{$reference} , #description
						\$reference_file, #reference alignment
						$ORF_RESULTS, #filehandle
						\%TCC, #path to config
						\%status_of_ORF_fasta, #contains all ORF fastas as keys
						\$crt_sample_filepath,
						\@suspicious_sequences
						);
					}
				}
			} else {
				my $crt_reference_fasta = File::Spec -> catfile ($TCC{'reference_fastas'}, $search_tool.'_AA_main_positive_company.fasta');
				print "\tjoining fastas to $crt_reference_fasta\n";
				unless (-s $crt_reference_fasta){ 
					&join_fastas(\%{$collection_contents{'main_positive'}{$search_tool}{'company'}}, \$crt_reference_fasta);
				}
				if ($search_tool eq 'jackhmmer'){
					#~ print "\t\t\t###running jackhmmer...";
					&run_jackhmmer(
						\'jackhmmer_AA_main_positive_company_collection', #name of the reference file
						\'jackhmmer_AA_main_positive_company', #description
						\$crt_reference_fasta , #reference alignment
						$ORF_RESULTS, #filehandle
						\%TCC, #path to config
						\%status_of_ORF_fasta, #contains all ORF fastas as keys
						\$crt_sample_filepath,
						\@suspicious_sequences
					);	
				}elsif ($search_tool eq 'blastp'){
					#~ print "\t\t\t###running blastp...";
					&run_blastp(
						\'blastp_AA_main_positive_company_collection', #name of the reference file
						\'blastp_AA_main_positive_company', #description
						\$crt_reference_fasta , #reference alignment
						$ORF_RESULTS, #filehandle
						\%TCC, #path to config
						\%status_of_ORF_fasta, #contains all ORF fastas as keys
						\$crt_sample_filepath,
						\@suspicious_sequences
					);	
				}elsif ($search_tool eq 'mmseqs'){
					#~ print "\t\t\t###running mmseqs...";
					&run_mmseqs(
						\'mmseqs_AA_main_positive_company_collection', #name of the reference file
						\'mmseqs_AA_main_positive_company', #description
						\$crt_reference_fasta , #reference alignment
						$ORF_RESULTS, #filehandle
						\%TCC, #path to config
						\%status_of_ORF_fasta, #contains all ORF fastas as keys
						\$crt_sample_filepath,
						\@suspicious_sequences
					);	
				}
			}
			
			
		}
	}
	###checking suspicious sequences vs NR via blastp
#	print Dumper @suspicious_sequences;
	print '#' x 70,"\n";
	my $crt_basename = basename($crt_sample);
	$crt_basename =~ s/\.[^\.]+$//;
	(my $crt_sample_ID = $crt_basename) =~ s/_ORFs_AA//;
	if (@suspicious_sequences){
		print "$crt_sample_ID: ", scalar @suspicious_sequences," sequence(s)!\n";
		my $crt_ORF_CSV = File::Spec -> catfile ($TCC{'ORF_dir'}, $crt_sample_ID.'_ORF_data.csv');
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
		my %suspicious_ORFs_AA; #collects the AA sequences of the ORFs of the sample
		my $crt_sample_suspicious_NT_fasta = File::Spec -> catfile ($TCC{'result_dir'}, $crt_sample_ID.'_suspicious_NT.fasta');
			open (my $CRT_SUSPICIOUS_NT, '>', $crt_sample_suspicious_NT_fasta) or die "could not write to '$crt_sample_suspicious_NT_fasta': $!\n";
			my $crt_sample_suspicious_fasta = File::Spec -> catfile ($TCC{'result_dir'}, $crt_sample_ID.'_suspicious_AA.fasta');
			open (my $CRT_SUSPICIOUS, '>', $crt_sample_suspicious_fasta) or die "could not write to '$crt_sample_suspicious_fasta': $!\n";
			my $crt_sample_suspicious_blasted = File::Spec -> catfile ($TCC{'result_dir'}, $crt_sample_ID.'_suspicious_AA_blasted.csv');
			print "writing to '$crt_sample_suspicious_fasta'...\n";
			foreach my $sequence_header (@suspicious_sequences){
				print "\t\t$sequence_header\n";
				my $crt_full_seq = $sequence_header.','.$ORF_data{'>'.$sequence_header};
				print $CRT_SUSPICIOUS_NT ">$sequence_header\n$ORF_data{'>'.$sequence_header}\n";
				foreach my $crt_ORF (sort keys %{$ORF_data_by_sequence{$sequence_header}}){
					my @cols = split (',', $ORF_data_by_sequence{$sequence_header}{$crt_ORF});
					print $CRT_SUSPICIOUS ">$crt_ORF\n$cols[7]\n";
					$suspicious_ORFs_AA{$crt_ORF} = $cols[7];
				}
			}
			close ($CRT_SUSPICIOUS);
			close ($CRT_SUSPICIOUS_NT);
			###blasting all sample sequences vs the whole nr database
			print "blasting suspicious ORFs from $crt_sample_ID vs $TCC{'blastp_db'}...\n";
			my $t1  = time();
			unless ( $TCC{'only_evaluate'}){
				if ($TCC{'blastp_db'} =~ m/-remote/){
					#system(qq($TCC{'blastp'} -query $crt_sample_suspicious_fasta -out $crt_sample_suspicious_blasted $TCC{'blastp_settings'} -db $TCC{'blastp_db'} -outfmt \"10 qseqid qlen qframe sframe length pident nident qcovhsp qstart qend sstart send qseq sseq evalue score bitscore sacc\")) and die "Fatal: blast failed for '$crt_sample_suspicious_fasta': $!\n";
					print "running remote blast!\n";
					if (-s $crt_sample_suspicious_blasted) {
						system (qq(rm $crt_sample_suspicious_blasted));
					}
					my %crt_suspicious;
					&fasta_2_hash(\$crt_sample_suspicious_fasta,\%crt_suspicious);
					for my $seq (keys %crt_suspicious){		
						my $crt_content = "$seq\n$crt_suspicious{$seq}\n";
						my $crt_fasta = 'tmp.fasta';
						my $crt_out = 'tmp_blasted.csv';
						&write_file(\$crt_fasta, \$crt_content);
						system(qq($TCC{'blastp'} -query $crt_fasta -out $crt_out $TCC{'blastp_settings'} -db $TCC{'blastp_db'} -outfmt \"10 qseqid qlen qframe sframe length pident nident qcovhsp qstart qend sstart send qseq sseq evalue score bitscore sacc\"));
						if (-s $crt_out){
							system(qq(cat $crt_out >> $crt_sample_suspicious_blasted));
							system (qq(rm $crt_out));
						}
						system (qq(rm $crt_fasta));
					}
				} else {	
				system(qq($TCC{'blastp'} -query $crt_sample_suspicious_fasta -out $crt_sample_suspicious_blasted $TCC{'blastp_settings'} -db $TCC{'blastp_db'} -num_threads $TCC{'nCPU'} -outfmt \"10 qseqid qlen qframe sframe length pident nident qcovhsp qstart qend sstart send qseq sseq evalue score bitscore sacc\")) and die "Fatal: blast failed for '$crt_sample_suspicious_fasta': $!\n";
				}
			}
			
			my $t2  = time();
			my $tdiff = $t2 - $t1;
			my $hitcounter;
			print "reading BLAST output '$crt_sample_suspicious_blasted'...\n";
			open (my $BLAST_OUT, '<', $crt_sample_suspicious_blasted) or die "could not read from '$crt_sample_suspicious_blasted': $!\n";
			while (my $line = <$BLAST_OUT>){
				chomp $line;
				next if ($line =~ m/^#/);
				next if ($line =~ m/^$/);
				my @cols = split (',', $line);
				$cols[-1] =~ s/\.\d\s*$//; #get rid of version number
				$hitcounter++;
				print $ORF_RESULTS "$crt_sample_ID,$cols[0],$cols[-1],nr,$tdiff,blastp_vs_nr,$cols[-4],$cols[-3],$cols[-2]\n";
			}
			if ($hitcounter){
				print "\t\tfound $hitcounter hit(s)!\n";
			} else {
				print "\t\tno hits found\n";
				print $ORF_RESULTS "$crt_sample_ID,NA,NA,nr,$tdiff,blastp_vs_nr,NA,NA,NA\n";
			}
	} else {
		print "no suspicious sequences in sample!\n";
		print $ORF_RESULTS "$crt_sample_ID,NA,NA,nr,NA,blastp_vs_nr,NA,NA,NA\n";
	}
		#~ $status_of_ORF_fasta{$crt_sample} 
	print $ORF_FASTAS_DONE "$crt_sample,$status_of_ORF_fasta{$crt_sample}\n";
		
	
}
	
	
print "TRAVIS Core finished successfully!\n";	
	

###subroutines
sub join_fastas {
	my $href_mother_hash = $_[0]; 
	my $sref_outfile = $_[1];
	my %sequences;
	foreach my $reference (keys %{$href_mother_hash}){
		foreach my $reference_file (@{$href_mother_hash -> {$reference}}){
			print "\t\t\t$reference_file\n";
			open (my $FASTA, '<', $reference_file) or print "could not read from '$reference_file': $!\n";
			my $current_header;
			my $skip_flag = 0;
			while (my $line = <$FASTA>){
				chomp $line;
				if ($line =~ m/^>/){
					#~ print "before: $line\n";
					my @header_elements= split ('__', $line);
					pop @header_elements unless ((scalar @header_elements <= 1) or ($TCC{'reference_library'} eq 'blind'));
					$current_header = join ('__', @header_elements);
					#~ print "after: $current_header\n";
					if ($sequences{$current_header}) {
						$skip_flag = 1;
					} else {
						$skip_flag = 0;
					}
				#~ print "header: $line\n";
				} else {
					$sequences{$current_header} .= $line unless ($skip_flag);
				}
			}	 
		}
	}
	open (my $FH, '>',$$sref_outfile ) or die "could not write to '$$sref_outfile : $!\n";
	foreach my $seq (keys %sequences){
		print  $FH "$seq\n$sequences{$seq}\n";
	}
	close ($FH);
 }
 
 
 
 
sub write_file {
	#takes a file path and a string to create a file with the string as content
	my $sref_file = $_[0]; #file path
	my $sref_content = $_[1]; #content to be written into the file
	open (my $FH, '>',$$sref_file) or die "could not write to '$$sref_file': $!\n";
	print $FH $$sref_content;
	close $FH;
	#~ print "$$sref_file written\n";
}


sub  slurp_fasta {
	my $sref_fasta_file = $_[0];
	my $sref_seq_container = $_[1];
	open (my $FASTA, '<', $$sref_fasta_file) or print "could not read from '$$sref_fasta_file': $!\n";
	my $current_header;
	while (my $line = <$FASTA>){
		$$sref_seq_container  .= $line;
	} 
}


sub fasta_2_hash {
	my $sref_fasta_file = $_[0];
	my $href_seq_container = $_[1];
	open (my $FASTA, '<', $$sref_fasta_file) or print "could not read from '$$sref_fasta_file': $!\n";
	my $current_header;
	while (my $line = <$FASTA>){
		chomp $line;
		if ($line =~ m/^>/){
			$current_header = $line;
			#~ print "header: $line\n";
		} else {
			$href_seq_container -> {$current_header} .= uc $line;
		}
	} 
	print "hashed $$sref_fasta_file...\n";
}

sub run_mmseqs{
	#~ return;
	my $sref_mmseqs_db_name = $_[0];
	my $sref_description = $_[1];
	my $sref_fasta = $_[2]; #file path
	my $sFHref_ORF_RESULTS = $_[3];
	my $href_TCC = $_[4];
	my $href_status_of_ORF_fasta = $_[5];
	my $sref_crt_sample = $_[6];
	my $aref_suspicious_sequences = $_[7];
	
	unless (-s $$sref_fasta ){
		print "$$sref_fasta does not exist or is empty...no need to run a search\n";
		return;
	}
	
	my $mmseqs_db_path = File::Spec -> catfile ($href_TCC -> {'reference_fastas'}, $$sref_mmseqs_db_name.'.mmseqsdb');
	my $mmseqs_results = File::Spec -> catdir ($href_TCC -> {'result_dir'}, 'mmseqs_results');
	&check_dir(\$mmseqs_results);
	my $mmseqs = $href_TCC -> {'mmseqs'};
	my $mmseqs_search_settings  = $href_TCC -> {'mmseqs_search_settings'};
	unless (-s $mmseqs_db_path){
		print "\tbuilding mmseqs database '$mmseqs_db_path'...\n";
		(system(qq($mmseqs createdb $$sref_fasta $mmseqs_db_path)) and die "Fatal: could not build mmseqs database from '$$sref_fasta'\: $!\n") unless ($href_TCC -> {'only_evaluate'});
	} else {
		print "\t'$mmseqs_db_path' already exists...\n";
	}
	#~ my $blastp = $href_TCC -> {'blastp'};
	my $ncpu = $href_TCC -> {'nCPU'};
	#~ my $blastp_settings = $href_TCC -> {'blastp_settings'};
	
	my $hitcounter;
	(my $mmseqs_query_db_path = $$sref_crt_sample) =~ s/fasta$/mmseqsdb/;
	unless (-s $mmseqs_query_db_path){
		print "\tbuilding mmseqs database '$mmseqs_query_db_path'...\n";
		my $mmseqs = $href_TCC -> {'mmseqs'};
		(system(qq($mmseqs createdb $$sref_crt_sample $mmseqs_query_db_path)) and die ("Fatal: could not build mmseqs database from '$$sref_crt_sample'\: $!\n")) unless ($href_TCC -> {'only_evaluate'});
	} else {
		print "\t'$mmseqs_query_db_path' already exists...\n";
	}
	my $crt_basename = basename($$sref_crt_sample);
	$crt_basename =~ s/\.[^\.]+$//;
	(my $crt_sample_ID = $crt_basename) =~ s/_ORFs_AA//;
	my $crt_tmp_results = File::Spec -> catdir ($mmseqs_results, $$sref_mmseqs_db_name.'_vs_'.$crt_basename.'_tmp');
	&check_dir(\$crt_tmp_results);
	my $mmseqs_result_db_path = File::Spec -> catdir ($mmseqs_results, $$sref_mmseqs_db_name.'_vs_'.$crt_basename.'.mmseqsdb');
	if (-s $mmseqs_result_db_path){
		(system(qq(rm $mmseqs_result_db_path))) unless ($href_TCC -> {'only_evaluate'});
	}
	#~ #print "$crt_basename\n";
	my $crt_result = File::Spec -> catfile ($mmseqs_results, $$sref_mmseqs_db_name.'_vs_'.$crt_basename.'.tsv');
	#~ #	print "$crt_result\n";
	my $t1  = time();
	unless ( $href_TCC -> {'only_evaluate'}){
		system(qq(rm -r $crt_tmp_results/*));
		system(qq($mmseqs search $mmseqs_query_db_path $mmseqs_db_path $mmseqs_result_db_path $crt_tmp_results --threads $ncpu $mmseqs_search_settings)) and die "Fatal: mmseqs failed for '$$sref_crt_sample': $!\n";
		print "search done...converting DBs...\n";
		system(qq($mmseqs convertalis $mmseqs_query_db_path $mmseqs_db_path $mmseqs_result_db_path $crt_result --threads $ncpu)) and die "Fatal: mmseqs conversion failed for '$$sref_crt_sample': $!\n";
	}
	my $t2 = time();
	my $tdiff = $t2 - $t1;

	open (my $CRT_RESULT, '<', $crt_result) or die "could not read from '$crt_result': $!\n";
	while (my $line = <$CRT_RESULT>){
		chomp $line;
		next if ($line =~ m/^#/);
		next if ($line =~ m/^$/);
		my @cols = split (/\t/, $line);
		$hitcounter++;
		print $sFHref_ORF_RESULTS "$crt_sample_ID,$cols[0],$cols[1],$$sref_mmseqs_db_name,$tdiff,mmseqs,$cols[-2],$cols[-1],NA\n";
		$cols[0] =~ m/(.+)_(ORF_\d+)/;
		my $seq_pt1 = $1;
		push(@$aref_suspicious_sequences, $seq_pt1) unless ($seq_pt1 ~~ @$aref_suspicious_sequences);
	}
	if ($hitcounter){
		print "\t\tfound $hitcounter hit(s)!\n";
		
		$href_status_of_ORF_fasta -> {basename($$sref_crt_sample)} = 'positive';
	} else {
		print "\t\tno hits found\n";
		print $sFHref_ORF_RESULTS "$crt_sample_ID,NA,NA,$$sref_mmseqs_db_name,$tdiff,mmseqs,NA,NA,NA\n";
	}
	
}




sub run_blastp{
	#~ return;
	my $sref_blast_db_name = $_[0];
	my $sref_description = $_[1];
	my $sref_fasta = $_[2]; #file path
	my $sFHref_ORF_RESULTS = $_[3];
	my $href_TCC = $_[4];
	my $href_status_of_ORF_fasta = $_[5];
	my $sref_crt_sample = $_[6];
	my $aref_suspicious_sequences = $_[7];
	
	unless (-s $$sref_fasta ){
		print "$$sref_fasta does not exist or is empty...no need to run a search\n";
		return;
	}
	
	my $blast_db_path = File::Spec -> catfile ($href_TCC -> {'reference_fastas'}, $$sref_blast_db_name.'.blastdb');
	my $blast_results = File::Spec -> catdir ($href_TCC -> {'result_dir'}, 'blast_results');
	&check_dir(\$blast_results);
	unless (-s $blast_db_path.'.phr'){
		print "\tbuilding blast database '$blast_db_path'...\n";
		my $makeblastdb = $href_TCC -> {'makeblastdb'};
		(system(qq($makeblastdb -in $$sref_fasta -dbtype prot -out $blast_db_path -title $blast_db_path )) and die "Fatal: could not build blast database from '$$sref_fasta '\: $!\n") unless ($href_TCC -> {'only_evaluate'});
	} else {
		print "\t'$blast_db_path' already exists...\n";
	}
	my $blastp = $href_TCC -> {'blastp'};
	my $ncpu = $href_TCC -> {'nCPU'};
	my $blastp_settings = '-evalue 10';
	$blastp_settings = $href_TCC -> {'blastp_settings'} if ($href_TCC -> {'blastp_settings'});
	
	my $hitcounter;
	my $crt_basename = basename($$sref_crt_sample);
	$crt_basename =~ s/\.[^\.]+$//;
	(my $crt_sample_ID = $crt_basename) =~ s/_ORFs_AA//;
	#~ #print "$crt_basename\n";
	my $crt_result = File::Spec -> catfile ($blast_results, $$sref_blast_db_name.'_vs_'.$crt_basename.'.csv');
	#	print "$crt_result\n";
	my $t1  = time();
	unless ( $href_TCC -> {'only_evaluate'}){
		(system(qq($blastp -query $$sref_crt_sample -out $crt_result $blastp_settings -db $blast_db_path -num_threads $ncpu -outfmt \"10 qseqid qlen qframe sframe length stitle pident nident qcovhsp qstart qend sstart send qseq sseq evalue score bitscore sacc\")) and die "Fatal: blast failed for '$$sref_crt_sample': $!\n")  unless ($href_TCC -> {'only_evaluate'});
	}
	my $t2  = time();
	my $tdiff = $t2 - $t1;

	open (my $CRT_RESULT, '<', $crt_result) or die "could not read from '$crt_result': $!\n";
	while (my $line = <$CRT_RESULT>){
		chomp $line;
		next if ($line =~ m/^#/);
		next if ($line =~ m/^$/);
		my @cols = split (',', $line);
		$cols[5] =~ s/\s+$//;
		$hitcounter++;
		print $sFHref_ORF_RESULTS "$crt_sample_ID,$cols[0],$cols[5],$$sref_blast_db_name,$tdiff,blastp,$cols[-4],$cols[-3],$cols[-2]\n";
		$cols[0] =~ m/(.+)_(ORF_\d+)/;
		my $seq_pt1 = $1;
		push(@$aref_suspicious_sequences, $seq_pt1) unless ($seq_pt1 ~~ @$aref_suspicious_sequences);
	}
	if ($hitcounter){
		print "\t\tfound $hitcounter hit(s)!\n";
		$href_status_of_ORF_fasta -> {basename($$sref_crt_sample)} = 'positive';
	} else {
		print "\t\tno hits found\n";
		print $sFHref_ORF_RESULTS "$crt_sample_ID,NA,NA,$$sref_blast_db_name,$tdiff,blastp,NA,NA,NA\n";
	}
	
}

sub run_jackhmmer{
	#~ return;
	my $sref_reference_name = $_[0];
	my $sref_description = $_[1];
	my $sref_fasta = $_[2]; #file path
	my $FH_ORF_RESULTS = $_[3];
	my $href_TCC = $_[4];
	my $href_status_of_ORF_fasta = $_[5];
	my $sref_crt_sample = $_[6];
	my $aref_suspicious_sequences = $_[7];
	
	unless (-s $$sref_fasta ){
		print "$$sref_fasta does not exist or is empty...no need to run a search\n";
		return;
	}
	
	my $jackhmmer_results = File::Spec -> catdir ($href_TCC -> {'result_dir'}, 'jackhmmer_results');
	&check_dir(\$jackhmmer_results);
	my $jackhmmer = $href_TCC -> {'jackhmmer'};
	my $ncpu = $href_TCC -> {'nCPU'};
	my $jackhmmer_settings = '-E 10';
	$jackhmmer_settings = $href_TCC -> {'jackhmmer_settings'} if ($href_TCC -> {'jackhmmer_settings'});
	
	my $hitcounter;
	my $crt_basename = basename($$sref_crt_sample);
	$crt_basename =~ s/\.[^\.]+$//;
	(my $crt_sample_ID = $crt_basename) =~ s/_ORFs_AA//;
	my $crt_result = File::Spec -> catfile ($jackhmmer_results, $$sref_reference_name.'_vs_'.$crt_basename);
	my $t1  = time();
	unless ( $href_TCC -> {'only_evaluate'}){
		(system(qq($jackhmmer $jackhmmer_settings --cpu $ncpu -o /dev/null --tblout $crt_result $$sref_fasta $$sref_crt_sample)) and die "Fatal: jackhmmerfailed for '$$sref_crt_sample': $!\n") unless ($href_TCC -> {'only_evaluate'}) ;
	}
	my $t2  = time();
	my $tdiff = $t2 - $t1;
	
	open (my $CRT_RESULT, '<', $crt_result) or die "could not read from '$crt_result': $!\n";
	while (my $line = <$CRT_RESULT>){
		chomp $line;
		next if ($line =~ m/^#/);
		$hitcounter++;
		$line =~ s/\h+/\t/g;
		my @cols = split (/\t/, $line);
			#~ print "$cols[0]\n";
		print $FH_ORF_RESULTS "$crt_sample_ID,$cols[0],$cols[2],$$sref_reference_name,$tdiff,jackhmmer,$cols[4],$cols[5],NA\n";
		$cols[0] =~ m/(.+)_(ORF_\d+)/;
		my $seq_pt1 = $1;
		push(@$aref_suspicious_sequences, $seq_pt1) unless ($seq_pt1 ~~ @$aref_suspicious_sequences);
	}
	if ($hitcounter){
		print "\t\tfound $hitcounter hit(s)!\n";
		$href_status_of_ORF_fasta -> {basename($$sref_crt_sample)} = 'positive';
	} else {
		print "\t\tno hits found\n";
		print $FH_ORF_RESULTS "$crt_sample_ID,NA,NA,$$sref_reference_name,$tdiff,jackhmmer,NA,NA,NA\n";
	}
}




sub run_hmmer{
	
	my $sref_hmm_name = $_[0];
	my $sref_description = $_[1];
	my $sref_alignment = $_[2]; #file path
	my $sFHref_ORF_RESULTS = $_[3];
	my $href_TCC = $_[4];
	my $href_status_of_ORF_fasta = $_[5];
	my $sref_crt_sample = $_[6];
	my $aref_suspicious_sequences = $_[7];
	
	unless (-s $$sref_alignment ){
		print "$$sref_alignment does not exist or is empty...no need to run a search\n";
		return;
	}
	
	my $hmm_path = File::Spec -> catfile ($href_TCC -> {'reference_fastas'}, $$sref_hmm_name.'.hmm');
	my $hmmsearch_results = File::Spec -> catdir ($href_TCC -> {'result_dir'}, 'hmmsearch_results');
	&check_dir(\$hmmsearch_results);
	unless (-s $hmm_path){
		print "\tbuilding hmm '$$sref_hmm_name'...\n";
		my $hmmbuild = $href_TCC -> {'hmmbuild'};
		(system(qq($hmmbuild --informat afa -n $$sref_hmm_name $hmm_path $$sref_alignment )) and die "Fatal: could not build hmm from '$$sref_alignment '\: $!\n") unless ($href_TCC -> {'only_evaluate'});
	} else {
		print "\t'$hmm_path' already exists...\n";
	}
	my $hmmsearch = $href_TCC -> {'hmmsearch'};
	my $ncpu = $href_TCC -> {'nCPU'};
	my $hmmsearch_settings = '-E 10';
	$hmmsearch_settings = $href_TCC -> {'hmmsearch_settings'} if ($href_TCC -> {'hmmsearch_settings'});
	my $hitcounter;
	my $crt_basename = basename($$sref_crt_sample);
	$crt_basename =~ s/\.[^\.]+$//;
	(my $crt_sample_ID = $crt_basename) =~ s/_ORFs_AA//;
	my $crt_result = File::Spec -> catfile ($hmmsearch_results, $$sref_hmm_name.'_vs_'.$crt_basename);
	my $t1  = time();
	unless ( $href_TCC -> {'only_evaluate'}){
		(system(qq($hmmsearch $hmmsearch_settings --cpu $ncpu -o /dev/null --tblout $crt_result $hmm_path $$sref_crt_sample)) and die "Fatal: hmmsearch failed for '$$sref_crt_sample': $!\n") unless ($href_TCC -> {'only_evaluate'});
	}
	my $t2  = time();
	my $tdiff = $t2 - $t1;
	open (my $CRT_RESULT, '<', $crt_result) or die "could not read from '$crt_result': $!\n";
	while (my $line = <$CRT_RESULT>){
		chomp $line;
		next if ($line =~ m/^#/);
		$line =~ s/\h+/\t/g;
		my @cols = split (/\t/, $line);
		$hitcounter++;
		#~ #print "$cols[0]\n";
		print $sFHref_ORF_RESULTS "$crt_sample_ID,$cols[0],$$sref_description,$$sref_hmm_name,$tdiff,hmmer,$cols[4],$cols[5],NA\n";
		$cols[0] =~ m/(.+)_(ORF_\d+)/;
		my $seq_pt1 = $1;
		push(@$aref_suspicious_sequences, $seq_pt1) unless ($seq_pt1 ~~ @$aref_suspicious_sequences);
	}
	if ($hitcounter){
		print "\t\tfound $hitcounter hit(s)!\n";
		$href_status_of_ORF_fasta -> {basename($$sref_crt_sample)} = 'positive';
	} else {
		print "\t\tno hits found\n";
		print $sFHref_ORF_RESULTS "$crt_sample_ID,NA,NA,$$sref_hmm_name,$tdiff,hmmer,NA,NA,NA\n";
	}
	
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


sub ORF_data_2_fasta{
	#converts the ORF_data table to a fasta that only contains the ORF ID and the AA sequence
	my $seq_header=0;
	my $seq_position=8;
	my $sref_table=$_[0];
	my $sref_fasta=$_[1];
	my $header='no';
	open (my $TABLE, '<', $$sref_table) or die "could not read from '$$sref_table': $!\n";
	open (my $FASTA, '>', $$sref_fasta) or die "could not write to '$$sref_fasta': $!\n";
	while (my $line = <$TABLE>){
		if ($header=~ m/yes/i){
			$header='no';
			next;
		}
		chomp $line;
		my @cols = split (',',$line);
		if ($cols[$seq_position]) {
			print $FASTA ">$cols[$seq_header]\n";
			print $FASTA $cols[$seq_position], "\n";
		}
	}
	print "converted $$sref_table to $$sref_fasta\n";
}


sub get_sample_ORF_data {
	#loops through each sequence of a fasta or a fastq or a fastq.gz
	my $sref_in = $_[0];
	my $sref_out = $_[1];
	my $sref_format = $_[2];
	my $sref_min_ORF_length = $_[3];
	my $sref_max_ORF_length = $_[4];
	if ($$sref_format eq 'fasta'){
		print "reading fasta: '$$sref_in'...\n";
		open (my $FH_IN, '<',$$sref_in) or die "could not read from '$$sref_in': $!\n";
		open (my $FH_OUT, '>',$$sref_out) or die "could not write to '$$sref_out': $!\n";
		my $crt_header;
		my $crt_seq;
		while (my $line = <$FH_IN>){
			chomp $line;
			if ($line =~ m/^>(.+)/){
				my $crt_seq_ID = $1;
				if (!$crt_header){
					($crt_header = $crt_seq_ID) =~ s/[^a-zA-Z0-9_\-]/_/g;
				} else {
					#~ print "evaluating $crt_header...\n";
					my $collected_ORF_data = &get_ORF_data(\$crt_header,\$crt_seq, $sref_min_ORF_length, $sref_max_ORF_length);
					if ($collected_ORF_data){
						print $FH_OUT ">$crt_header,$crt_seq\n";
						print $FH_OUT $collected_ORF_data;
					} 
					($crt_header = $crt_seq_ID) =~ s/[^a-zA-Z0-9_\-]/_/g;
					undef $crt_seq;
				}
			} else {
				$crt_seq .= $line;
			}
		}
		my $collected_ORF_data = &get_ORF_data(\$crt_header,\$crt_seq, $sref_min_ORF_length, $sref_max_ORF_length);
		if ($collected_ORF_data){
			print $FH_OUT ">$crt_header,$crt_seq\n";
			print $FH_OUT $collected_ORF_data;
		} 
		close ($FH_IN);
		close ($FH_OUT);
	} elsif ($$sref_format eq 'fasta.gz'){
		print "reading fasta.gz is not yet implemented :/ '$$sref_in'...\n";
	} elsif ($$sref_format eq 'fastq'){
		print "reading fastq: '$$sref_in'...\n";
		open (my $FH_IN, '<',$$sref_in) or die "could not read from '$$sref_in': $!\n";
		open (my $FH_OUT, '>',$$sref_out) or die "could not write to '$$sref_out': $!\n";
		my $crt_header;
		my $crt_seq;
		my $seq_count;
		my $linecount;
		my %identicals;
		my $number_of_different_sequences;
		while (my $line = <$FH_IN>){
			chomp $line;
			$linecount++;
			#~ last if ($linecount >= 40000);
			
			if ($linecount % 4 == 1){
				$line =~ m/^@(.+)/;
				($crt_header = $1) =~ s/[^a-zA-Z0-9_\-]/_/g;
				$seq_count++;
				if ($seq_count % 10000 == 0){
					print "evaluated $seq_count sequences\n";
				}
				#~ print "header: $crt_header \n";
			} elsif ($linecount % 4 == 2){
				#~ print "seq; $line\n";
				if (!$identicals{$line}){
					$number_of_different_sequences++;
					my $collected_ORF_data = &get_ORF_data(\$crt_header,\$line, $sref_min_ORF_length, $sref_max_ORF_length);
					if ($collected_ORF_data){
						print $FH_OUT ">$crt_header,$line\n";
						print $FH_OUT $collected_ORF_data;
					}
				} 
				push (@{$identicals{$line}}, $crt_header);
			}
		}
		my $collected_ORF_data = &get_ORF_data(\$crt_header,\$crt_seq, $sref_min_ORF_length, $sref_max_ORF_length);
		if ($collected_ORF_data){
			print $FH_OUT ">$crt_header,$crt_seq\n";
			print $FH_OUT $collected_ORF_data;
		} 
		close ($FH_IN);
		close ($FH_OUT);
		
		(my $identicals = $$sref_out) =~ s/ORF_data\.csv/identical_sequences.csv/;
		open ($FH_OUT, '>',$identicals) or die "could not write to '$identicals': $!\n";	
		foreach my $key (keys %identicals){
			my $entries = join('&', @{$identicals{$key}});
			my $number_of_entries = scalar @{$identicals{$key}};
			print $FH_OUT "$number_of_entries,$key,$entries\n";
		}
			
		close ($FH_OUT);	
		print "$number_of_different_sequences different sequences in total of $seq_count sequences\n";
			
			
			
	} elsif ($$sref_format eq 'fastq.gz'){
		print "reading fastq.gz is not yet implemented :/ '$$sref_in'...\n";
	} else {
		die "'$$sref_format' is not supportet!\n";
	}
	print "wrote $$sref_out\n";
}


sub get_ORF_data{
	#takes references to a sequence ID, a DNA string and a minimum ORF length (AA).
	#looks for ORFs in all reading frames that have lengths above or equal to the minimum ORF length
	#returns a reference to a comma separated value table containing ORF_ID, reading frame, reading frame orientation, start position on the original DNA sequence,
	# end position on the original DNA sequence, length of the ORF in NT, length of the ORF in AA, the NT sequence of the ORF, the AA sequence of the ORF
	my $sref_sequence_ID = $_[0];
	my $sref_original_sequence = $_[1];
	my $sref_minimal_ORF_length = $_[2];
	my $sref_maximal_ORF_length = $_[3];
	print "a minimal ORF length below 10 AA will probably cause troubles!\n" if ($$sref_minimal_ORF_length < 10);
	my $ORF_data;
	my $number_of_ORFs;
	#~ my @reading_frame_orientations = ('forward', 'reverse_complement');
	for (my $reading_frame = 1; $reading_frame <= 3; $reading_frame++){ #loop through each reading frame (1->3)
		my $sequence = uc $$sref_original_sequence;
		foreach my $reading_frame_orientation ('forward', 'reverse_complement'){ #for each direction
			#~ print "processing $reading_frame_orientation frame $reading_frame of $sequence...\n";
			my $original_sequence_length = length($sequence);
			if ($reading_frame_orientation eq 'reverse_complement'){
				$sequence = reverse($sequence);
				$sequence =~ tr/ABCDGHMNRSTUVWXY/TVGHCDKNYSAABWXR/;
			}
			$sequence = substr($sequence, $reading_frame - 1);
			my $translated_sequence;
			&translate(\$sequence, \$translated_sequence);
			my @ORFs = split('\*',$translated_sequence);
			foreach my $ORF (@ORFs){
				my $ORF_AA_length = length($ORF);
				next if ($ORF_AA_length < $$sref_minimal_ORF_length);
				next if ($ORF_AA_length > $$sref_maximal_ORF_length);
				$number_of_ORFs++;
				my $crt_ORF_ID =join ('_', "$$sref_sequence_ID\_ORF", sprintf ("%03d", $number_of_ORFs));
				my $ORF_AA_offset = index($translated_sequence, $ORF);
				my $ORF_NT = substr($sequence, 3 * $ORF_AA_offset , 3 * $ORF_AA_length);
				my $ORF_NT_start = 3 * $ORF_AA_offset + $reading_frame;
				$ORF_NT_start = $original_sequence_length - $ORF_NT_start + 1 if ($reading_frame_orientation eq 'reverse_complement');
				my $ORF_NT_length = length($ORF_NT);
				my $ORF_NT_end = $ORF_NT_start + $ORF_NT_length - 1 ;
				$ORF_NT_end = $ORF_NT_start - $ORF_NT_length + 1 if ($reading_frame_orientation eq 'reverse_complement');
				$ORF_data .= "$crt_ORF_ID,$reading_frame_orientation,$reading_frame,$ORF_NT_start,$ORF_NT_end,$ORF_NT_length,$ORF_AA_length,$ORF_NT,$ORF\n";
			}
		}
	}
	return $ORF_data;
}


sub translate{
	#takes a DNA string and retuns its translation
	my $sref_DNA = $_[0];
	my $sref_translated_sequence = $_[1]; #in undef, out def
	my %codon_translations = (
	    'TAA' => '*',
	    'TAG' => '*',
	    'TGA' => '*',
	    'TCC' => 'S', 
	    'TTC' => 'F', 
	    'TCA' => 'S',
	    'TCT' => 'S', 
	    'TCG' => 'S',
	    'TGC' => 'C',
	    'TTA' => 'L',
	    'TAC' => 'Y',
	    'TTG' => 'L',
	    'TTT' => 'F',
	    'TAT' => 'Y',
	    'TGG' => 'W',
	    'CTA' => 'L',
	    'CTC' => 'L',
	    'TGT' => 'C', 
	    'CTG' => 'L',
	    'CTT' => 'L', 
	    'CCA' => 'P',
	    'CCG' => 'P',
	    'CCT' => 'P',
	    'CAC' => 'H',
	    'CCC' => 'P',
	    'CAA' => 'Q',
	    'CGA' => 'R',
	    'CAT' => 'H',
	    'CGT' => 'R',
	    'CGC' => 'R',
	    'ATA' => 'I',
	    'CAG' => 'Q',
	    'ATC' => 'I',
	    'ATT' => 'I',
	    'CGG' => 'R',
	    'ATG' => 'M',
	    'ACA' => 'T',
	    'ACT' => 'T',
	    'ACG' => 'T',
	    'AAT' => 'N',
	    'AAA' => 'K',
	    'AGG' => 'R',
	    'AAC' => 'N',
	    'ACC' => 'T',
	    'AGC' => 'S',
	    'AAG' => 'K',
	    'AGA' => 'R',
	    'AGT' => 'S',
	    'GTA' => 'V',
	    'GTC' => 'V',
	    'GTT' => 'V',
	    'GCC' => 'A',
	    'GCA' => 'A',
	    'GTG' => 'V',
	    'GCT' => 'A',
	    'GAT' => 'D',
	    'GCG' => 'A',
	    'GAC' => 'D',
	    'GGG' => 'G',
	    'GAG' => 'E',
	    'GGA' => 'G',
	    'GGT' => 'G',
	    'GAA' => 'E',
	    'GGC' => 'G'
	);
	for (my $i =0; $i <= (length($$sref_DNA) -3); $i += 3){
		my $crt_codon = substr($$sref_DNA, $i,3);
		my $crt_translation = 'X';
		$crt_translation =$codon_translations{$crt_codon} if ($codon_translations{$crt_codon}); 
		$$sref_translated_sequence .= $crt_translation;
	}
}
