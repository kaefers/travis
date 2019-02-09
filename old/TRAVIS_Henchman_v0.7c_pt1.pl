#!/usr/bin/env perl
use warnings;
use strict;
use LWP::Simple;
use autodie;
use File::Path qw(make_path);	
use File::Basename;
use File::Spec;
use Data::Dumper;
use Net::FTP;

#~ unless (@ARGV){
	#~ die "usage: perl TRAVIS_Henchman_v0.6c.pl <TRAVIS_control_center.csv>\n";
#~ }
#~ my $TCC_file = shift @ARGV;
my $TCC_file = '/home/skaefer/sync/work/reoviridae/reo_mini_henchman_splitter/reo_mini_test_TCC.csv';
#~ my $TCC_file = '/home/skaefer/sync/work/parkinson/synuclein_TCC.csv';

print '#' x 70,"\n";
print '#' x 70,"\n";
print "\t\t\tThis is TRAVIS Henchman v0.7c pt1\n";
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


###load local reference database into a hash
#local reference database is a CSV with ID,description,sequence,origin,name,acronym,genus
#1st dimension: column name starting with description, 2nd dimension: ID
my %local_reference_database_entries;
if (-s $TCC{'local_reference_database'}){
	open (my $LOC_REF_DB, '<', $TCC{'local_reference_database'}) or die "could not read from '$TCC{'local_reference_database'}': $!\n";
	print "reading $TCC{'local_reference_database'}...\n";
	while (my $line = <$LOC_REF_DB>){
		chomp $line;
		next if ($line =~ m/fail/); #skip entries that failed previously
		my @cols = split (',', $line);
		$local_reference_database_entries{'descriptions'}{$cols[0]} = $cols[1];
		$local_reference_database_entries{'sequences'}{$cols[0]} = $cols[2];
		$local_reference_database_entries{'origins'}{$cols[0]} = $cols[3];
		$local_reference_database_entries{'virus_names'}{$cols[0]} = $cols[4];
		$local_reference_database_entries{'virus_acronyms'}{$cols[0]} = $cols[5];
		$local_reference_database_entries{'virus_genera'}{$cols[0]} = $cols[6];
	}
	close ($LOC_REF_DB);
}
###

#~ my $look_at_this;


###read the reference library and (down)load the data
my $all_NT_col; #will contain the index of the column called 'All_NT_ACC'
my $main_NT_col; #will contain the index of the column called '<MAIN>_NT_ACC'
my $main_PID_col; #will contain the index of the column called '<MAIN>_PID'
my $all_PID_col; #will contain the index of the column called 'All_PID'
my $main_name; #will contain thename of the main protein/gene/whatever '<MAIN>_NT_ACC'
my @header_elements; #will contain the names of the columns of the reference library
my @header_cols; #will contain the indices of columns that will become the sequence header in the specified order in $TCC{'header_names'}
my @split_by_cols; #will contain the indices of columns by which the main references will be split specified by $TCC{'split_references'}
my @failed; #will collect failed entries for later review
my %sequences; #multidimensional hash that will contain sequences. 1st dimension: NT/AA, 2nd dimension: main/company, 3rd dimension: split groups


###open and loop through each line of the reference library
#download and sort all data
#create a preparsed library for faster re-runs
unless ($TCC{'use_preparsed'}){ #set preparsed to 'yes' if it has not been declared in TCC. 'yes' is the default
	$TCC{'use_preparsed'} = 'yes';
}
(my $preparsed_reference_library = $TCC{'reference_library'}) =~ s/\.csv$/_preparsed.csv/;
unless (-s $preparsed_reference_library){
	print "$preparsed_reference_library is either not existing or empty! I will create one!\n";
	$TCC{'use_preparsed'} = 'no';
}
	
my $REFLIB; #File handle to be used for reading the reference library. either the original or the preparsed version depending on $TCC{'use_preparsed'}
my $PREPARSED_REFLIB; #File handle to be used for writing the preparsed reference library
if ($TCC{'use_preparsed'} eq 'no'){
	print "will read from $TCC{'reference_library'} \n";
	open ($REFLIB, '<', $TCC{'reference_library'}) or die "could not read from '$TCC{'reference_library'}': $!\n";
	open ($PREPARSED_REFLIB, '>', $preparsed_reference_library) or die "could not write to '$preparsed_reference_library': $!\n";
} else {
	open ($REFLIB, '<', $preparsed_reference_library) or die "could not read from '$preparsed_reference_library': $!\n";
	print "will read from $preparsed_reference_library \n";
	
}

&check_dir(\$TCC{'reference_fastas'}); #to see if this directory exists an files can be written within

my $c;
while (my $line = <$REFLIB>){
	chomp $line;
	#~ print "$line\n";
	$c++;
	#~ last if ($c >2);
	my @cols = split (',', $line);
	###check reference library table headers
	unless (@header_elements){ #do this if @header_elements is not yet defined
		@header_elements = @cols;
		###get the accession number column numbers by looping through them
		#example: All_NT_ACC,RdRp_NT_ACC,RdRp_PID
		for (my $i=0; $i < scalar @cols; $i++){
			if ($cols[$i] eq 'all_NT_ACC'){
				$all_NT_col= $i;
			} elsif ($cols[$i] =~ m/(.*)_NT_ACC/){
				$main_NT_col= $i;
				$main_name = $1;
			} elsif (($main_NT_col) and ($cols[$i] =~ m/$main_name\_PID/)){
				$main_PID_col= $i;
				$all_PID_col = $main_PID_col + 1;
			}
		}
		###get order of fasta sequence header elements
		foreach my $put_in (split ('&', $TCC{'header_names'})){
			for (my $i=0; $i < scalar @cols; $i++){
				if ($cols[$i] eq $put_in){
					push (@header_cols, $i);
				}
			}
		}
		###
		#get columns of split factors
		if ($TCC{'split_references'}){
			foreach my $split_by  (split('&', $TCC{'split_references'})){
				for (my $i=0; $i < scalar @cols; $i++){
					if ($cols[$i] eq $split_by){
						push (@split_by_cols, $i);
					}
				}
			}
		}
		###
		if ($PREPARSED_REFLIB) { #put this header into the preparsed reference library if it is to be created
			print $PREPARSED_REFLIB $line, ",all_PID\n";
		}
		next;
	}
	###
	#~ print '#' x 70,"\n";
	print "retrieving data for $cols[1]...\n";
	###create header according to $TCC{'header_names'}
	my $crt_header_elements;
	foreach my $colnumber (@header_cols){
		$crt_header_elements .= '__'.$cols[$colnumber];
	}
	#print "$crt_header_elements\n";
	###
	###load all accessions and corresponding annotated proteins
	my $crt_all_NT_ACCs = $cols[$all_NT_col]; #should contain either a list of IDs separated by '&' or a link to an NCBI assembly report
	my @collected_NT_ACCs; #will contain all retrieved NT accessions for the current reference library line to be put into the preparsed reference library
	my @collected_PIDs; #will contain all retrieved PIDs for the current reference library line to be put into the preparsed reference library
	(my $main_NT_ACC = $cols[$main_NT_col]) =~ s/\.\d\s*$//; #single NT accession and get rid of version number
	(my $main_PID  = $cols[$main_PID_col]) =~ s/\.\d\s*$//; #single PID that can be found in $main_NT_ACC and get rid of version number
	my @crt_all_NT_ACCs; #$crt_all_NT_ACCs parsed into an array of accession numbers:
	if ($crt_all_NT_ACCs =~ m/ftp.ncbi.nlm.nih.gov/){
		for (my $tries = 1; $tries <=3; $tries++){
			last if (@crt_all_NT_ACCs);
			print "downloading assembly report from NCBI (try $tries of 3)...\n";
			&read_assembly_report(
				\$crt_all_NT_ACCs, #ftp path
				\@crt_all_NT_ACCs #undefined array
			);
		}
		unless (@crt_all_NT_ACCs){
			push(@failed, $crt_all_NT_ACCs);
			print "fuckit on crt_all_NT_ACCs: '$crt_all_NT_ACCs' in 177!\n";
		}
	} elsif ($crt_all_NT_ACCs eq 'NA'){
		print "no other sequences? sure?\n";
	} else {
		print "parsing provided accession numbers...\n";
		@crt_all_NT_ACCs = split ('&', $crt_all_NT_ACCs);		
	}
	push (@collected_NT_ACCs, @crt_all_NT_ACCs); #actually collect
	###loop through all NT accession numers, get associated PIDs, download sequences and additional info
	###will be stored in %sequences and %local_reference_database_entries
	my $crt_sequence; #will contain a sequence
	my $crt_description; #will contain NCBI description of a sequence
	print "loading data...\n";
	my $main_or_company; #a switch for 'Main' and 'Company' depending on the current NT_ACC or PID
	foreach my $crt_NT_ACC (@crt_all_NT_ACCs){
		#~ print "next NT_ACC $crt_NT_ACC\n";
		$crt_NT_ACC =~ s/\.\d$//; #get rid of version number
		if ($crt_NT_ACC eq $main_NT_ACC){
			$main_or_company =  $main_name;
		} else {
			$main_or_company = 'company';
		}
		
		#~ print "for: $crt_NT_ACC\n";
		my @crt_PIDs; #will contain all PIDs associated with $crt_NT_ACC
		if ($cols[$all_PID_col]){ #if there is content in $cols[$all_PID_col] (i.e. if we are reading a preparsed reference library):
			@crt_PIDs = split ('&', $cols[$all_PID_col]); #take the PIDs from that variable
		} else { 
			&fetch_PIDs( #get it from the genebank entry
				\$crt_NT_ACC, #NT accession number
				\@crt_PIDs #undefined array to be filled with PIDs
			);
		}
		
		###see if the sequence for $crt_NT_ACC is already in the local database. 
		#if not, retrieve it from NCBI
		if (!$local_reference_database_entries{'sequences'}{$crt_NT_ACC}){
			&fetch_sequence(
				\$crt_NT_ACC, #NT accession number
				\$crt_description, #NCBI description of the sequence to be retrieved
				\$crt_sequence, #sequence to be retrieved
				\'nuccore' #which database
			);
			if (!$crt_sequence) { #check if efetch worked. if it failed:
				push(@failed, $crt_NT_ACC); #store the $crt_NT_ACC in @failed for later review
				print "fuckit on crt_sequence: '$crt_sequence' in 220!\n";
				next; #head to next sequence
			}
			###store data in %local_reference_database_entries
			$local_reference_database_entries{'descriptions'}{$crt_NT_ACC}= $crt_description; #freshly downloaded
			$local_reference_database_entries{'sequences'}{$crt_NT_ACC}= $crt_sequence; #freshly downloaded
			$local_reference_database_entries{'origins'}{$crt_NT_ACC}=  'self'; #it's a NT-sequence. so it has its origin in itself
			$local_reference_database_entries{'virus_names'}{$crt_NT_ACC}= $cols[1]; #from reference library
			$local_reference_database_entries{'virus_acronyms'}{$crt_NT_ACC}= $cols[0]; #from reference library
			$local_reference_database_entries{'virus_genera'}{$crt_NT_ACC}= $cols[4]; #from reference library
			###
		} else {
			$crt_sequence  = $local_reference_database_entries{'sequences'}{$crt_NT_ACC};
			$crt_description = $local_reference_database_entries{'descriptions'}{$crt_NT_ACC};
		}
		
		if ($main_or_company ne 'company'){
			print "$crt_NT_ACC: NT says main\n";
			$sequences{'NT'}{$main_or_company}{'all'} .= ">$crt_NT_ACC$crt_header_elements\n$crt_sequence\n"; #sort the entry into the 'all' group...
		} else {
			$sequences{'NT'}{$main_or_company}{'all'} .= ">$crt_NT_ACC$crt_header_elements\__$crt_description\n$crt_sequence\n"; #sort the entry into the 'all' group...
		}
		#..and into all split groups:
		foreach my $split_by_col (@split_by_cols){
			$sequences{'NT'}{$main_or_company}{$header_elements[$split_by_col].'_'.$cols[$split_by_col]} .= ">$crt_NT_ACC$crt_header_elements\__$crt_description\n$crt_sequence\n";
		}
		if (@crt_PIDs) {
			push (@collected_PIDs, @crt_PIDs) ; #actually collect 
		} else {
			print "no PID found on $crt_description\n";
		}
		###
		###set some values of the local reference database for the associated PIDs
		foreach my $crt_PID(@crt_PIDs){
			$crt_PID =~ s/\.\d\s*$//; #get rid of version number
			$local_reference_database_entries{'origins'}{$crt_PID}=  $crt_NT_ACC; #it's a PID found in $crt_NT_ACC
			$local_reference_database_entries{'virus_names'}{$crt_PID}= $cols[1]; #from reference library
			$local_reference_database_entries{'virus_acronyms'}{$crt_PID}= $cols[0]; #from reference library
			$local_reference_database_entries{'virus_genera'}{$crt_PID}= $cols[4]; #from reference library
		}
		###
	}
	###get the sequences for all associated PIDs
	next unless (@collected_PIDs); 	
	my @PIDs_done;
	foreach my $crt_PID (@collected_PIDs){	
		$crt_PID =~ s/\.\d\s*$//; #get rid of version number
		next if ($crt_PID ~~ @PIDs_done);
		next if ($crt_PID =~ m/^\s*$/);
		push (@PIDs_done, $crt_PID);
		#~ print "crt_PID: $crt_PID\n";
		if ($crt_PID eq $main_PID){
			$main_or_company = $main_name;
		} else {
			$main_or_company = 'company';
		}
		if (!$local_reference_database_entries{'sequences'}{$crt_PID}){
			&fetch_sequence(
				\$crt_PID, #PID
				\$crt_description, #NCBI description of the sequence to be retrieved
				\$crt_sequence, #sequence to be retrieved
				\'protein' #which database
			);
			if (!$crt_sequence) { #check if efetch worked. if it failed:
				print "fuckit on crt_sequence: '$crt_sequence'  by crt_PID $crt_PID in 281!\n";
				push(@failed, $crt_PID); #store the $crt_PID in @failed for later review
				next; #head to next sequence
			}
			
			###store data in %local_reference_database_entries
			$local_reference_database_entries{'descriptions'}{$crt_PID}= $crt_description; #freshly downloaded
			$local_reference_database_entries{'sequences'}{$crt_PID}= $crt_sequence; #freshly downloaded
		###
		} else {
			$crt_sequence  = $local_reference_database_entries{'sequences'}{$crt_PID};
			$crt_description = $local_reference_database_entries{'descriptions'}{$crt_PID};
		}	
		if ($main_or_company ne 'company'){
			print "$crt_PID is Main! $crt_description\n";
			$sequences{'AA'}{$main_or_company}{'all'} .= ">$crt_PID$crt_header_elements\__$main_name\n$crt_sequence\n";#sort the entry into the 'all' group...
		} else {
			#~ print  "$crt_description\n";
			$sequences{'AA'}{$main_or_company}{'all'} .= ">$crt_PID$crt_header_elements\__$crt_description\n$crt_sequence\n";#sort the entry into the 'all' group...
		}
		#..and into all split groups:
		foreach my $split_by_col (@split_by_cols){
			$sequences{'AA'}{$main_or_company}{$header_elements[$split_by_col].'_'.$cols[$split_by_col]} .= ">$crt_PID$crt_header_elements\__$main_name\n$crt_sequence\n";
		}
	}
	###
	###if the preparsed reference library is to be written
	if ($PREPARSED_REFLIB){
		for (my $i = 0; $i <= scalar @cols; $i++){ #loop through all columns
			if ($i == $all_NT_col){ #identify the all NT column
				print $PREPARSED_REFLIB join ('&', @collected_NT_ACCs), ','; #print all NT Accs joined by '&' in the column
			} elsif ($i == $all_PID_col){ #identifz the all PID column
				print $PREPARSED_REFLIB join ('&', @collected_PIDs), "\n"; #print all PIDs joined by '&' in the column
			} else {
				print $PREPARSED_REFLIB $cols[$i], ','; #otherwise print the original column content
			}
		}
	}
	###
	
	
	
	###check if the main accessions ere really part of all provided accessions
	##later......
	#unless ($local_reference_database_entries{'sequences'}{$main_PID} and $local_reference_database_entries{'sequences'}{$main_NT_ACC}){
		#~ my @alligot = (@all_crt_NT_IDs, @all_crt_AA_IDs);
		#~ $lookatthis .= "###################\n$cols[1]\n";
		#~ foreach my $whatigot (@alligot){
			#~ $lookatthis .= "$whatigot: $descriptions{$whatigot} in $origins{$whatigot}\n";
		#~ }
		#~ print "could not find the specified main NT or AA sequence. Please review\n";
		#~ print "######################\n";
		#~ next;
	#~ }
	###
}
###

###update local database
###later: make path dirname of $TCC{'local_reference_database'}
open (my $LOC_REF_DB, '>', $TCC{'local_reference_database'}) or die "can not write to '$TCC{'local_reference_database'}': $!\n";
foreach my $ID (keys %{$local_reference_database_entries{'sequences'}}){
	print $LOC_REF_DB "$ID,$local_reference_database_entries{'descriptions'}{$ID},$local_reference_database_entries{'sequences'}{$ID},$local_reference_database_entries{'origins'}{$ID},$local_reference_database_entries{'virus_names'}{$ID},$local_reference_database_entries{'virus_acronyms'}{$ID},$local_reference_database_entries{'virus_genera'}{$ID}\n";
}
print "updated local database '$TCC{'local_reference_database'}'\n";
###
###create all group fastas and calculate alignments of main sequences

my %alignments_of_fastas; #contains the path to an aligned version of the fasta. 'none' if the fasta has not been aligned because there is only one sequence inside
foreach my $seq_type (keys %sequences){
	foreach my $set_type (keys %{$sequences{$seq_type}}){
		foreach my $group (keys %{$sequences{$seq_type}{$set_type}}){
			print "$seq_type $set_type $group\n";
			print " $set_type $group\n";
			my $crt_fasta =  File::Spec -> catfile($TCC{'reference_fastas'}, $TCC{'database_name'}."_$seq_type\_$set_type\_$group.fasta");
			
			#~ &add_to_header(
			$sequences{$seq_type}{$set_type}{$group} =~ s/(>[^\n]+)\n/$1\___$set_type\_$group\n/g;
			#~ print $sequences{$seq_type}{$set_type}{$group}. "\n";	
			#~ )
			
			print "writing $crt_fasta\n";
			&write_file(
				\$crt_fasta, #file path to fasta
				\$sequences{$seq_type}{$set_type}{$group} #sequence(s) in fasta format
			);
			#~ print "$sequences{$seq_type}{$set_type}{$group}\n\n\n";
			if ($set_type ne 'company'){

				(my $crt_aln_fasta = $crt_fasta) =~ s/\.fasta$/_aln.fasta/;
				my $crt_fasta_basename = basename($crt_fasta);	
				my $crt_aln_fasta_basename = basename($crt_aln_fasta);	
				$alignments_of_fastas{$crt_fasta} = $crt_aln_fasta;
				 
				my $how_many_sequences = `grep '>' $crt_fasta | wc -l` ;
				chomp $how_many_sequences;
				print "$how_many_sequences sequences!\n";
				if ($how_many_sequences > 1){
					
					if ($seq_type eq 'AA') {
						#~ print "aligning...\n";
						$TTT_content .= "$set_type\_$group,main\_$main_name,$crt_fasta_basename,$crt_aln_fasta_basename,$how_many_sequences,all,hmmer&jackhmmer&mmseqs&blastp,on\n"; #add  an entry to TTT
					
					
						#~ if ($TCC{'mafft_binaries'}){
							#~ system(qq(export MAFFT_BINARIES="$TCC{'mafft_binaries'}";$TCC{'mafft'} $TCC{'mafft_settings'} --thread $TCC{'nCPU'} $crt_fasta > $crt_aln_fasta)) ;
						#~ }else {
							#~ system(qq($TCC{'mafft'} $TCC{'mafft_settings'} --thread $TCC{'nCPU'} $crt_fasta > $crt_aln_fasta)) ;
						#~ }
						
					}
				} else {
					if ($seq_type eq 'AA') {
						$TTT_content .= "$set_type\_$group,main\_$main_name,$crt_fasta_basename,NA,$how_many_sequences,all,jackhmmer&mmseqs&blastp,on\n"; #add  an entry to TTT
					}
					#system(qq(cat $crt_fasta > $crt_aln_fasta)) ;
				}
			} 
			else {
				$alignments_of_fastas{$crt_fasta} = 'none';
			}
		}
	}
}
###
&write_file(\$TCC{'TTT'},\$TTT_content); #write TTT

if (@failed){
	print "failed to retrieve ", scalar @failed	, " items:\n";
	foreach my $fail (@failed){
		print $fail, "\n";
	}
}


print "TRAVIS Henchman pt1 completed\n";


exit;













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



if (@failed){
	print "failed to retrieve ", scalar @failed	, " items:\n";
	foreach my $fail (@failed){
		print $fail, "\n";
	}
}


print "TRAVIS Henchman pt1 completed\n";

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



sub fetch_sequence{
	#download NCBI efetch data in fasta format for a given NT accession number and parse it for description and sequence
	my $sref_accession_number = $_[0]; #accession number or PID
	my $sref_crt_description = $_[1]; #in: defined or undefined, out: defined
	my $sref_crt_sequence = $_[2]; #in: defined or undefined, out: defined
	my $sref_database = $_[3]; #'nuccore' or 'protein'
	my $gb_entry = get "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$$sref_database&id=$$sref_accession_number&rettype=fasta";
	if (!$gb_entry){
		print "failed to retrieve sequence for $$sref_accession_number\n";
	} else {
		undef $$sref_crt_sequence;
		my @lines = split("\n", $gb_entry);
		for (my $l=1; $l < scalar @lines;$l++){
			chomp $lines[$l];
			$$sref_crt_sequence .= $lines[$l];
		}
		($$sref_crt_description = $lines[0]) =~ s/[^a-zA-Z0-9 _\-\[\]]|$$sref_accession_number[^\s]* |\s*\[[^\]]+\]//g;
		$$sref_crt_description =~ s/\s/_/g;
	}
	sleep(1);  #NCBI does not like too many failed attempts in a short amount of time ;)
}




sub fetch_PIDs {
	#download NCBI efetch data in genebank XML format for a given NT accession number and parse it for associated PIDs
	#PIDs will be stored in an array
	my $sref_accession_number=$_[0]; #a NT accession number
	my $aref_crt_PIDs = $_[1]; #in: undefined, out: defined
	print "fetching PIDs for $$sref_accession_number from NCBI...\n";
	#http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001542&rettype=gb&retmode=xml
	my $gb_entry = get "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$$sref_accession_number&rettype=gb&retmode=xml"; #downloads the entry into a string
	if (!$gb_entry){
		sleep(1); 
		print "failed to retrieve genebank entry for $$sref_accession_number!\n";
		#~ push(@$aref_crt_PIDs, 'failed');
	} else {
		my @gb_entry_lines= split("\n",$gb_entry); #split the entry to read it line by line
		for (my $i=1; $i < scalar @gb_entry_lines; $i++){
			#<GBQualifier_name>protein_id</GBQualifier_name>
			# <GBQualifier_value>NP_056796.1</GBQualifier_value>
			if ($gb_entry_lines[$i] =~ m/<GBQualifier_name>protein_id<\/GBQualifier_name>/){ #look for the protein_id-tag
				$i++; #skip to next line
				#~ print "PID found!\n";
				$gb_entry_lines[$i] =~ m/<GBQualifier_value>([^<]+)<\/GBQualifier_value>/; #extract the PID
				(my $crt_PID= $1) =~ s/\.\d\s*$//; #get rid of version number
				push(@$aref_crt_PIDs, $crt_PID) unless ($crt_PID =~ m/^\s*$/); #storesPID
			}
		}
		sleep(1); 
	}
}

sub read_assembly_report {
	#downloads an assembly report from the NCBI ftp-server and parses it for accession numbers. 
	my $sref_assembly_report = $_[0]; #ftp path
	my $aref_crt_all_NT_ACCs = $_[1]; #in: undefined, out: defined
	$$sref_assembly_report =~ s/ftp\:\/\/ftp.ncbi.nlm.nih.gov//; #extract file path from ftp path
	my $ftp =  new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1) or die "Cannot connect to ftp.ncbi.nlm.nih.gov: $!"; #connect to ftp-server
	$ftp -> login; $ftp->binary; #login to ftp-server
	$ftp -> get ($$sref_assembly_report) or print "get failed ", $ftp->message; #download file
	$ftp -> quit; #disconnect
	$$sref_assembly_report = basename($$sref_assembly_report);	
	###read the assembly report. Interesting is the tab-separated table after the commented lines
	open (my $ASSEMBLY_REPORT, '<', $$sref_assembly_report) or return "could not read from '$$sref_assembly_report': $!\n";
	while (my $line = <$ASSEMBLY_REPORT>){
		chomp $line;
		next if ($line =~ m/^#/); #skip comments
		my @cols = split("\t", $line); #split line by tab
		#~ print $cols[6], "\n";
		push (@$aref_crt_all_NT_ACCs , $cols[6]); #$cols[6] contains the accession numbers, pushed to the array
	}
	close ($ASSEMBLY_REPORT);
	system(qq(rm $$sref_assembly_report)); #remove file 
	sleep(1); 
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


