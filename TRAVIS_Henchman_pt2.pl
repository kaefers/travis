#!/usr/bin/env perl
use warnings;
use strict;
use File::Path qw(make_path);	
use File::Basename;
use File::Spec;
use Data::Dumper;


unless (@ARGV){
	die "usage: perl TRAVIS_Henchman_pt2.pl <TRAVIS_control_center.csv>\n";
}
my $TCC_file = shift @ARGV;


print '#' x 70,"\n";
print '#' x 70,"\n";
print "\t\t\tThis is TRAVIS Henchman pt2 v20190209\n";
print '#' x 70,"\n";
print '#' x 70,"\n";
print "\n";
###reading TRAVIS Control Center
my %TCC; #contains configuration parameters from $TCC_file 
&read_TCC(
	\$TCC_file, #file path
	\%TCC #undef hash
);
print "connected to TRAVIS Control Center\n";
###
my $TTT_content;


print "reading $TCC{'TTT'}...\n";
open (my $TTT, '<', $TCC{'TTT'}) or die "could not read from '$TCC{'TTT'}': $!\n";
while (my $line = <$TTT>){
	chomp $line;
	my @cols = split(',', $line);
	print "$line\n";
	
	
	

						#######$TTT_content .= "$set_type\_$group,main\_$main_name,$crt_fasta_basename,NA,$how_many_sequences,all,jackhmmer&mmseqs&blastp,on\n"; #add  an entry to TTT
					
	my $crt_fasta = File::Spec -> catfile ($TCC{'reference_fastas'},$cols[2]);
	my $crt_aln_fasta = File::Spec -> catfile ($TCC{'reference_fastas'}, $cols[3]);
	$TTT_content .= $line."\n";
	next if ($crt_aln_fasta eq 'NA');
	print "aligning $crt_fasta...\n";
	if ($TCC{'mafft_binaries'}){
		system(qq(export MAFFT_BINARIES="$TCC{'mafft_binaries'}";$TCC{'mafft'} $TCC{'mafft_settings'} --thread $TCC{'nCPU'} $crt_fasta > $crt_aln_fasta)) ;
	}else {
		system(qq($TCC{'mafft'} $TCC{'mafft_settings'} --thread $TCC{'nCPU'} $crt_fasta > $crt_aln_fasta)) ;
	}


}

close($TTT);




###cluster all company  AA sequences using mmseqs2
my $fasta =  File::Spec -> catfile($TCC{'reference_fastas'}, $TCC{'database_name'}.'_AA_company_all.fasta'); #name of the AA fasta that contains all company sequences
print "clustering all company AA based on $fasta...\n";
#filenames based on fasta:
(my $crt_DB = $fasta) =~ s/fasta$/db/;
(my $crt_clustered = $fasta) =~ s/fasta$/clstr/;
(my $crt_clu = $fasta) =~ s/fasta$/clu/;
(my $crt_DB_tempdir = $fasta) =~ s/\.fasta$/_tmpdir\//;
(my $crt_clu_seq = $fasta) =~ s/fasta$/clu_seq/;

#~ system(qq(rm -r $crt_DB_tempdir /home/skaefer/sync/work/reoviridae/references/*clu*));
&check_dir(\$crt_DB_tempdir);
###mmseqs standard clustering pipe
system(qq($TCC{'mmseqs'} createdb $fasta $crt_DB));
system(qq($TCC{'mmseqs'} cluster $crt_DB $crt_clu $crt_DB_tempdir --threads $TCC{'nCPU'} $TCC{'mmseqs_cluster_settings'}));
system(qq($TCC{'mmseqs'} createseqfiledb $crt_DB $crt_clu $crt_clu_seq));
system(qq($TCC{'mmseqs'} result2flat $crt_DB $crt_DB $crt_clu_seq $crt_clustered));
###

open (my $CLUSTERED, '<', $crt_clustered) or die "could not read from '$crt_clustered': $!\n";
my $nextseq = 0; #flag for indicating next cluster OR new sequence
my $crt_cluster_name; #stores the current cluster name
my $crt_header; #stores the current sequence header
my %clusters; #stores sequences: key = $crt_cluster_name, value = concatenated fasta 
my %cluster_descriptions; #stores all information about the protein that is contained in the header by $crt_cluster_name
my %cluster_sizes; #stores the number of sequences per cluster

###read the clustered file:
while (my $line = <$CLUSTERED>){
	chomp $line;
	#~ print "line: $line\n";
	if (($nextseq == 0) and ($line =~ m/^>/)){
		#~ print "\t\tjust another seq in the cluster?\n";
		$crt_header = $line;
		$crt_header =~ s/^\s+|\s+$//; #get rid of spaces at the beginning and the end. just in case
		$nextseq = 1;
	} elsif ($nextseq ==1 and $line =~ m/^>(.*)/){
		$crt_cluster_name = $1;
		#~ print "\t\t################\n\t\tno, it's the next cluster: $crt_cluster_name!\n";
	} else {
		###if it actually is the header of another sequence within the cluster
		#~ print "\t\tnext sequence!: $crt_header:$line\n";
		$cluster_sizes{$crt_cluster_name}++; #+1 sequence for this cluster
		$clusters{$crt_cluster_name} .= "$crt_header\n$line\n"; #add the sequence to the cluster hash
		my @descriptions = split('__', $crt_header); #split the header for extracting protein description
		my $crt_description = $descriptions[-2]; #protein description is the second last entry in the header....at least it should be...
		$crt_description =~ s/^\s+|\s+$//; #get rid of spaces at the beginning and the end. just in case
		#~ print "desc: $crt_description\n";
		push (@{$cluster_descriptions{$crt_cluster_name}}, $crt_description) unless ($crt_description ~~ @{$cluster_descriptions{$crt_cluster_name}}); #collect descriptions of the single proteins. skip, if another protein had the same description
		$nextseq = 0;
	}
}
close ($CLUSTERED);
###

my $number_of_clusters; #count tne number of clusters. also used for naming them
my $max_cluster_size=0; #f.y.i: what the biggest cluster's size is	
my $small_clusters; #f.y.i : how many clusters were below the set minimal cluster size
my $unclustered; #a fasta-format string that collects sequences that are in a cluster <= minimal cluster siye

###loop through the clusters
foreach my $cluster (keys %clusters){
	my $crt_description = join('__', @{$cluster_descriptions{$cluster}}); #join the collected descriptions for the cluster to a string
	$max_cluster_size = $cluster_sizes{$cluster} if ($cluster_sizes{$cluster} > $max_cluster_size); #set maximum cluster size to current cluster size, if it is bigger than a known one
	if ($cluster_sizes{$cluster} <= $TCC{'minimal_cluster_size'}){ #if the cluster size is below minimal cluster size:
		$small_clusters++; #+1 to the small clusters counter
		$unclustered .= $clusters{$cluster} ; #add the sequence to the unclustered sequence collection
	} else { #if the cluster size is above minimal cluster size 
		$number_of_clusters++; #+1 to the cluster counter	
		#~ print "$cluster_sizes{$cluster} sequences\n";
		my $crt_out = File::Spec -> catfile( $TCC{'reference_fastas'}, $TCC{'database_name'}.'_AA_company_cluster_'.$number_of_clusters.'.fasta'); #create a path to a fasta for the current cluster
		#~ print "writing to $crt_out\n";
		#~ $sequences{$seq_type}{$set_type}{$group} =~ s/(>[^\n]+)\n/$1\___$set_type\_$group\n/g;
		#~ print "$clusters{$cluster}\n#######\n";
		
		$clusters{$cluster} =~ s/___company_all//g;
		$clusters{$cluster} =~ s/(>[^\n]+)\s*\n/$1\___company_cluster_$number_of_clusters\n/g;
		#~ print "fucking \$1 $1\n";
		#~ print "$clusters{$cluster}\n#######";
		&write_file(
			\$crt_out,
			\$clusters{$cluster}
		); #write the fasta
		
		(my $crt_out_aln = $crt_out) =~ s/\.fasta$/_aln.fasta/;; #create a path to the aligned fasta for the current cluster
		#and align it:
		my $crt_out_basename = basename($crt_out);
		my $crt_out_aln_basename = basename($crt_out_aln);
		
		if ($TCC{'mafft_binaries'}){
			system(qq(export MAFFT_BINARIES="$TCC{'mafft_binaries'}";$TCC{'mafft'} $TCC{'mafft_settings'} --thread $TCC{'nCPU'} $crt_out > $crt_out_aln)) ;
		}else {
			system(qq($TCC{'mafft'} $TCC{'mafft_settings'} --thread $TCC{'nCPU'} $crt_out > $crt_out_aln)) ;
		}
		
		$TTT_content .= 'company_cluster_'.$number_of_clusters.",company__$crt_description,$crt_out_basename,$crt_out_aln_basename,$cluster_sizes{$cluster},$TCC{'sample_subset'},hmmer&blastp&mmseqs&jackhmmer,on\n"; #add  an entry to TTT
	}
	#~ print "#" x 70;
	#~ print "\ncluster $cluster:\n$cluster_sizes{$cluster} sequence(s)\n$crt_description\n$clusters{$cluster}";
}
###
###create the unclustered reference fasta
if ($unclustered){
	my $crt_unclustered_out = File::Spec -> catfile ($TCC{'reference_fastas'} , $TCC{'database_name'}.'_AA_company_all_unclustered.fasta');
	$unclustered =~ s/(>[^\n]+)\n/$1\_unclustered\n/g;
	&write_file(\$crt_unclustered_out,\$unclustered);
	my $crt_unclustered_out_basename = basename($crt_unclustered_out);
	my $how_many_sequences = `grep '>' $crt_unclustered_out | wc -l`;
	chomp $how_many_sequences;
	$TTT_content .= "company_unclustered,company_unclustered,$crt_unclustered_out_basename,NA,$how_many_sequences,$TCC{'sample_subset'},blastp&mmseqs&jackhmmer,on\n"; #add  an entry to TTT
	#~ print "writing to $crt_unclustered_out\n";
}
print "$number_of_clusters cluster(s)\nmax size: $max_cluster_size sequences\n";
print "$small_clusters cluster(s) with less or equal than $TCC{'minimal_cluster_size'} sequence(s)\n" if ($small_clusters);
&write_file(\$TCC{'TTT'},\$TTT_content); #write TTT



###



print "re-running TRAVIS Henchman might change reference fasta composition and names! deleting all previously calculated reference databases...\n";
my $reference_databases_dir = File::Spec  -> catdir ($TCC{'result_dir'});
system (qq(rm $reference_databases_dir/*.blastdb* $reference_databases_dir/*.hmm $reference_databases_dir/*.mmseqsdb* $reference_databases_dir/*.blastdb*)) ;




print "TRAVIS Henchman pt2 completed\n";

###subroutines

sub write_file {
	#takes a file path and a string to create a file with the string as content
	my $sref_file = $_[0]; #file path
	my $sref_content = $_[1]; #content to be written into the file
	open (my $FH, '>',$$sref_file) or die "could not write to '$$sref_file': $!\n";
	print $FH $$sref_content;
	close $FH;
	#~ print "$$sref_file written\n";
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

sub  slurp_file{
	my $sref_file = $_[0];
	my $sref_container = $_[1];
	open (my $FILE, '<', $$sref_file) or print "could not read from '$$sref_file': $!\n";
	while (my $line = <$FILE>){
		$$sref_container  .= $line;
	} 
}
