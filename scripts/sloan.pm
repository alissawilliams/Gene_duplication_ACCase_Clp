###############################################################
#sloan.pm 
#A module of subroutines including some taken from Tisdall 2001
###############################################################

################################################################
#BEGIN open_file
# From example10-5.pl

# open_file
#
#   - given filename, set filehandle

sub open_file {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh;

    unless(open($fh, $filename)) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh;
}
#END open_file
###############################################################

################################################################
#BEGIN open_output

# open_output
#
#   - given filename, set filehandle

sub open_output {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh_output;

    unless(open($fh_output, ">$filename")) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh_output;
}
#END open_output
###############################################################

################################################################
#BEGIN open_output_append

#   - given filename, set filehandle

sub open_output_append {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh_output_append;

    unless(open($fh_output_append, ">>$filename")) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh_output_append;
}
#END open_output_append
###############################################################



###############################################################
# BEGIN file_to_array
#
# A subroutine to get data from a file given its filename

sub file_to_array {
	use strict;
	use warnings;

    my($filename) = @_;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}
#END file_to_array
###############################################################

###############################################################
# BEGIN file_to_string
#
# A subroutine to get data from a file given its filename

sub file_to_string {
	use strict;
	use warnings;

    my($filename) = @_;

    # Initialize variables
    my $filedata;

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }
    
    while (<GET_FILE_DATA>){
    	$filedata .= $_;
    }
    
    close GET_FILE_DATA;

    return $filedata;
}
#END file_to_string
###############################################################


###############################################################
#BEGIN arrays2hash
#a subroutine to take two arrays of equal size and convert them to a hash
#keys taken from first array, values from the second.
#Note that arrays must be passed by reference (\@array1, \@array2)

sub arrays2hash {
	use strict;
	use warnings;

	(my $keyarray, my $valuearray) = @_;
	if (scalar(@$keyarray) != scalar(@$valuearray)) {
		die "Arrays differ in size: Mismatched number of keys and values"; 
	}
	
	my %newhash = ( );
	
	@newhash{ @$keyarray } = @$valuearray;


	return (%newhash);

}

#END arrays2hash
###############################################################

###############################################################################
#get_fasta_names_and_seqs
#returns matching arrays of fasta heders and seqs given a fastname filename

sub get_fasta_names_and_seqs {
	use strict;
	use warnings;

	my ($inputfilename) = @_;
	my @fasta_names = ();
	my @fasta_seqs= ();

		   
	unless ( open(FILEDATA, $inputfilename) ) {
		print STDERR "Cannot open file \"$inputfilename\"\n\n"; #print error message
		exit; #exit the program
	}	

	my @filedata = <FILEDATA>; #Read the lines of the file into an array
	close FILEDATA;
	
	my $seq_count = 0; #this will be used to keep track of the number of sequences
	foreach my $line (@filedata){
		if ($line =~ /^>/) { #if the line is a header line (begins with ">")...
			if ($line =~ /^>.*[\w]+/){
				my $partialLine = substr ($&, 1);
				push (@fasta_names, $partialLine); #add that line to an array of fasta names
				push (@fasta_seqs, ""); #and add a new blank element to an array of sequences
				++$seq_count; #also increment our counter which keeps track of sequence number
			}
		}else { #if the line's not blank or a header, add it to the current sequence 
			$fasta_seqs[$seq_count-1] .= $line;
		}
	}
	for (my $i = 0; $i < scalar (@fasta_seqs); ++$i){
		$fasta_seqs[$i] =~s/\s//g;
	}
	
	return (\@fasta_names, \@fasta_seqs);

}

###############################################################################
###############################################################################
#get_fasta_names_and_seqs_old
#returns matching arrays of fasta heders and seqs given a fastname filename

sub get_fasta_names_and_seqs_old {
	use strict;
	use warnings;

	my ($inputfilename) = @_;
	my @fasta_names = ();
	my @fasta_seqs= ();

		   
	unless ( open(FILEDATA, $inputfilename) ) {
		print STDERR "Cannot open file \"$inputfilename\"\n\n"; #print error message
		exit; #exit the program
	}	

	my @filedata = <FILEDATA>; #Read the lines of the file into an array
	close FILEDATA;
	
	my $seq_count = 0; #this will be used to keep track of the number of sequences
	foreach my $line (@filedata){
		chomp $line;
		if ($line =~ /^\s*$/) {next;} #ignore line if it is blank
		
		elsif ($line =~ /^>/) { #if the line is a header line (begins with ">")...
			if ($line =~ /^>.*[\w]+/){
				my $partialLine = substr ($&, 1);
				push (@fasta_names, $partialLine); #add that line to an array of fasta names
				push (@fasta_seqs, ""); #and add a new blank element to an array of sequences
				++$seq_count; #also increment our counter which keeps track of sequence number
			}
		}	
		
		else { #if the line's not blank or a header, add it to the current sequence 
			$fasta_seqs[$seq_count-1] .= $line;
		}
		
		$fasta_seqs[$seq_count-1] =~s/\s//g; #remove all whitespace from the current  sequence
	}
	
	return (\@fasta_names, \@fasta_seqs);

}

###############################################################################

###############################################################
# BEGIN revcom 
#
# A subroutine to compute the reverse complement of DNA sequence

sub revcom {
	use strict;
	use warnings;

    my($dna) = @_;

    # First reverse the sequence
    my $revcom = reverse $dna;

    # Next, complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C, etc.
    $revcom =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;

    return $revcom;
}
#END revcom
###############################################################


###############################################################
# BEGIN codon2aa
#
# A subroutine to translate a DNA codon to an amino acid using hash lookup

sub codon2aa {
	use strict;
	use warnings;

    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    '---' => '-',    # gap
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
            return "?";
    }
}

#END codon2aa
###############################################################

###############################################################
# BEGIN dna2peptide 
#
# A subroutine to translate DNA sequence into a peptide

sub dna2peptide {
	use strict;
	use warnings;

    my($dna) = @_;

    use strict;
    use warnings;
    
    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($dna,$i,3) );
    }

    return $protein;
}
#END dna2peptide
###############################################################

###############################################################
# BEGIN gff_to_hash_of_sparse_arrays 
#a function to take a gff file and return a hash, with genome/sequence names as keys 
#and sparse arrays as values. The position of each sparse array is designated 1 if it is covered by 
#a feature in the gff file and is undefined otherwise.

#NOTE: INDEXING IS DONE STARTING AT "1" not "0". BE CAREFUL WITH INTERFACING WITH OTHER SUBROUTINES

sub gff_to_hash_of_sparse_arrays {
	use strict;
	use warnings;
	
	my $file = shift (@_) or die ("ERROR in gff_to_hash_of_sparse_array. No filename provided\n\n");
	my @fileData = file_to_array ($file);
	my %HoA; #hash of arrays to store an array for each named sequence in the gff file
	
	foreach my $line (@fileData){
		my @splitLine = split (/\t/, $line);
		for (my $i = $splitLine[3]; $i <= $splitLine[4]; ++$i){
			$HoA{$splitLine[0]}[$i] = 1;
		}
	}
	return %HoA;
}

###############################################################
#oneLetter2threeLetter_aa_code
#uses hash to convert single letter amino acid code to three letter code:

sub oneLetter2threeLetter_aa_code{
	my $oneLetter = shift @_ or die ("ERROR in oneLetter2threeLetter_aa_code\n: No amino acid code provided\n\n");
	length($oneLetter)==1 or die ("ERROR in oneLetter2threeLetter_aa_code\n: Specified amino acid code is longer than one letter\n\n");
    my(%aaCodeHash) = (
    
    'A' => 'Ala',    
    'C' => 'Cys',    
    'D' => 'Asp',    
    'E' => 'Glu',    
    'F' => 'Phe',    
    'G' => 'Gly',    
    'H' => 'His',    
    'I' => 'Ile',    
    'K' => 'Lys',    
    'L' => 'Leu',    
    'M' => 'Met',    
    'N' => 'Asn',    
    'P' => 'Pro',    
    'Q' => 'Gln',    
    'R' => 'Arg',    
    'S' => 'Ser',    
    'T' => 'Thr',    
    'V' => 'Val',    
    'W' => 'Trp',    
    'Y' => 'Tyr'
	);
	
	if (exists $aaCodeHash{uc($oneLetter)}){
		return $aaCodeHash{uc($oneLetter)};
	}else {
		return "???";
	}
	
}

###############################################################
###############################################################
# BEGIN gff_add_to_hash_of_sparse_arrays 
#takes a gff file and a reference to an existing hash of arrays.
#adds 1 to any position covered by a feature in the gff file

sub gff_add_to_hash_of_sparse_arrays {
	use strict;
	use warnings;
	my $file = shift (@_) or die ("ERROR in gff_add_to_hash_of_sparse_array. No filename provided\n\n");
	my @fileData = file_to_array ($file);
	my $HoAref = shift (@_) or die ("ERROR in gff_add_to_hash_of_sparse_array. No hash reference provided\n\n");
	my %HoA = %$HoAref;
	
	foreach my $line (@fileData){
		my @splitLine = split (/\t/, $line);
		for (my $i = $splitLine[3]; $i <= $splitLine[4]; ++$i){
			$HoA{$splitLine[0]}[$i] = 1;
		}
	}
	return %HoA;
}
###############################################################

###############################################################
# BEGIN gff_subtract_from_hash_of_sparse_arrays 
#takes a gff file and a reference to an existing hash of arrays.
#adds 1 to any position covered by a feature in the gff file

sub gff_subtract_from_hash_of_sparse_arrays {
	use strict;
	use warnings;
	my $file = shift (@_) or die ("ERROR in gff_subtract_from_hash_of_sparse_array. No filename provided\n\n");
	my @fileData = file_to_array ($file);
	my $HoAref = shift (@_) or die ("ERROR in gff_subtract_from_hash_of_sparse_array. No hash reference provided\n\n");
	my %HoA = %$HoAref;
	
	foreach my $line (@fileData){
		my @splitLine = split (/\t/, $line);
		if (exists $HoA{$splitLine[0]}){
			for (my $i = $splitLine[3]; $i <= $splitLine[4]; ++$i){
				$HoA{$splitLine[0]}[$i] = undef;
			}
		}
	}
	return %HoA;
}
###############################################################


###############################################################
#Begin gff_to_fasta
#takes a gff file and a corresponding fasta file and extracts 
#a seq for each gff line and returns it as multi-fasta formated string

sub gff_to_fasta{
	use strict;
	use warnings;
	
	my $gffFile = shift(@_) or die ("ERROR: gff_to_fasta. Insufficent arguments\n");
	my $inputFasta = shift(@_) or die ("ERROR: gff_to_fasta. Insufficent arguments\n");
	my @fileData = file_to_array($gffFile);
	my %HoAoA; #hash of array of arrays. Keys are seq names. Values are arrays of 2-element arrays (start and stop positions) 
	foreach my $line (@fileData){
		if ($line =~ /^\s*$/ || $line =~ /^#/){next;}
		my @splitLine = split (/\t/, $line);
		my @startEndArray;
		if ($splitLine[6] eq '-'){
			@startEndArray = ($splitLine[4],$splitLine[3]);
		}else{
			@startEndArray = ($splitLine[3],$splitLine[4]);
		}
		my $startEndArrayRef = \@startEndArray;
		push (@{$HoAoA{$splitLine[0]}}, $startEndArrayRef);
	}
	return fasta_seq_grab_multiple_with_header($inputFasta, \%HoAoA);
}

#end gff_to_fasta
###############################################################



###############################################################
# BEGIN coords_from_sparse_array 
#
# a subroutine to take a sparse array (containing some undefined positions), and return coordinates of all defined positions

sub coords_from_sparse_array {
	use strict;
	use warnings;

    my @inputArray = @_;
    my @outputArray;
   
    
    for(my $i=0; $i < scalar @inputArray ; ++$i) {
        if (defined $inputArray[$i]){
        	push (@outputArray, $i);
        }
    }

    return @outputArray;
}
#END coords_from_sparse_array
###############################################################

###############################################################
# BEGIN coords2ranges
#
# a subroutine to take an array of coordinates (passed as a reference) and return an array of 2 elements arrays,
#containing start and end positions for consecutive runs of coordinates.
# runs will be connected if gaps are < maxgap. Set maxgap to 0 to only consider coninuous runs.
#maxgap = 0 is also the default value if no maxgap is provided as an argument 
#Pass input array as a reference (1st argument) and $maxgap as a scalar (2nd argument)

sub coords2ranges {
	use strict;
	use warnings;

	my $coordsRef = shift @_; #array passed as reference
	my $maxgap;
	unless ($maxgap = shift @_){
		$maxgap = 0; 
	}
	my @coords = sort {$a <=> $b} @$coordsRef;
	my $range_cnt = 0;
	my (@ranges, $i);
	if (scalar(@coords) == 0){
		return @ranges;
	}else{
		$ranges[0][0] = $coords[0];
		for ($i = 1; $i < @coords; ++$i) {
		  if ($coords[$i] - $coords[$i - 1] > $maxgap + 1) {
		    $ranges[$range_cnt][1] = $coords[$i-1];
		    ++$range_cnt;
		    $ranges[$range_cnt][0] = $coords[$i];
		  }
		}
		$ranges[$range_cnt][1] = $coords[scalar(@coords)-1];
		return @ranges; #an array of 2-element arrays
	}
}
#END coords2ranges
###############################################################
################################################################################

#get_file_names
#given directory, return an array of all file names
sub get_file_names {
	use strict;
	use warnings;

    my ($directory) = @_;
    my @files = (  );
    my @filedata =(  );
    	

    # Open the directory
    unless(opendir(DIRECTORY, $directory)) {
        print "Cannot open directory $directory!\n";
        exit;
    }
    
    # Read the directory, ignoring special entries starting with "."
    @files = grep (!/^\./, readdir(DIRECTORY));
    
    closedir(DIRECTORY);
    
    return (@files);
   
}
###############################################################################


###############################################################################
#fasta_seq_grab
#given filename, start position, end position, and (optionally) seqName, returns specified sequence

sub fasta_seq_grab {
	use strict;
	use warnings;

	my $usage = "\n\nERROR in fasta_seq_grab: arguments (inputFastaFile, StartPosition, EndPosition, SequenceName[optional:assumed to be first sequence if ommitted])\n\n";

	my $inputFile = shift(@_) or die ($usage);
	my $startPos = shift(@_) or die ($usage);
	my $endPos = shift(@_) or die ($usage);
	my $seqNameProvided = 0;
	my $seqName;
	
	if ($_[0]){
		$seqName = $_[0];
		$seqNameProvided = 1;
	}

	my ($seqNamesRef, $seqSeqsRef) = get_fasta_names_and_seqs($inputFile);
	my @seqNames = @$seqNamesRef;
	my @seqSeqs = @$seqSeqsRef;
	my %seqFastaHash =  arrays2hash (\@seqNames, \@seqSeqs);
	
	unless ($seqNameProvided){
		$seqName = $seqNames[0];
	}
	
	my $extractedSeq;


	if (exists $seqFastaHash{$seqName}){
		if ($endPos >= $startPos){
			$extractedSeq = substr($seqFastaHash{$seqName}, $startPos - 1, $endPos - $startPos + 1);	
		}elsif ($startPos > $endPos){
			$extractedSeq = revcom (substr($seqFastaHash{$seqName}, $endPos - 1, $startPos - $endPos + 1));
		}else{
			die ("\n\nERROR: Coordinates do not appear to be numeric\n\n");
		}
	}else {
		die ("\n\nERROR: The following sequence name \"$seqName\" was not found in the fasta file: $inputFile\n\n");
	}
	
	return $extractedSeq;

}
###############################################################################

###############################################################################
#fasta_seq_grab_with_header
#given filename, start position, end position, and (optionally) seqName, returns specified sequence
#Same as fasta_seg_grab except that it also returns fasta header in the form:
#  >SeqName_StartPosition-EndPosition

sub fasta_seq_grab_with_header {
	use strict;
	use warnings;

	my $usage = "\n\nERROR in fasta_seq_grab_with_header: arguments (inputFastaFile, StartPosition, EndPosition, SequenceName[optional:assumed to be first sequence if ommitted])\n\n";

	my $inputFile = shift(@_) or die ($usage);
	my $startPos = shift(@_) or die ($usage);
	my $endPos = shift(@_) or die ($usage);
	my $seqNameProvided = 0;
	my $seqName;
	
	if ($_[0]){
		$seqName = $_[0];
		$seqNameProvided = 1;
	}

	my ($seqNamesRef, $seqSeqsRef) = get_fasta_names_and_seqs($inputFile);
	my @seqNames = @$seqNamesRef;
	my @seqSeqs = @$seqSeqsRef;
	my %seqFastaHash =  arrays2hash (\@seqNames, \@seqSeqs);
	
	unless ($seqNameProvided){
		$seqName = $seqNames[0];
	}
	
	my $extractedSeq;


	if (exists $seqFastaHash{$seqName}){
		if ($endPos >= $startPos){
			$extractedSeq = substr($seqFastaHash{$seqName}, $startPos - 1, $endPos - $startPos + 1);	
		}elsif ($startPos > $endPos){
			$extractedSeq = revcom (substr($seqFastaHash{$seqName}, $endPos - 1, $startPos - $endPos + 1));
		}else{
			die ("\n\nERROR: Coordinates do not appear to be numeric\n\n");
		}
	}else {
		die ("\n\nERROR: The following sequence name \"$seqName\" was not found in the fasta file: $inputFile\n\n");
	}
	
	my $fastaFormat = ">$seqName"."_$startPos"."-$endPos\n"."$extractedSeq\n";
	
	return $fastaFormat;

}
###############################################################################

###############################################################################
#fasta_seq_grab_multiple_with_header
#similar to fasta_seq_grab. But takes a hash of array of arrays (e.g. as produced within gff_to_fasta) 
#key = seqname (corresponding to fasta headers). Each value is an array of 2 element arrays with start and stop positions
#makes mulitple simulataneous extractions from a fasta file and returns a fasta formated string
#output fasta header in the form:
#  >SeqName_StartPosition-EndPosition

sub fasta_seq_grab_multiple_with_header {
	use strict;
	use warnings;

	my $usage = "\n\nERROR in fasta_seq_grab_with_header: arguments (inputFastaFile, HashReference )\n\n";

	my $inputFile = shift(@_) or die ($usage);
	my $HoAoAref = shift(@_) or die ($usage);
	my %HoAoA = %$HoAoAref;
		
	my %seqFastaHash =  arrays2hash (get_fasta_names_and_seqs($inputFile));
	
	my $extractedSeq;
	my $outputString;

	foreach my $seqName (sort keys %HoAoA){
		my @AoA = @{$HoAoA{$seqName}};
		foreach my $arrayRef (@AoA){
			my ($startPos, $endPos) = @$arrayRef;
			if (exists $seqFastaHash{$seqName}){
				if ($endPos >= $startPos){
					$extractedSeq = substr($seqFastaHash{$seqName}, $startPos - 1, $endPos - $startPos + 1);	
				}elsif ($startPos > $endPos){
					$extractedSeq = revcom (substr($seqFastaHash{$seqName}, $endPos - 1, $startPos - $endPos + 1));
				}else{
					die ("\n\nERROR: Coordinates do not appear to be numeric\n\n");
				}
				$outputString .= ">$seqName"."_$startPos"."-$endPos\n"."$extractedSeq\n";
			}else {
				die ("\n\nERROR: The following sequence name \"$seqName\" was not found in the fasta file: $inputFile\n\n");
			}
		}
	}
	
	return $outputString;

}
###############################################################################


###############################################################################
#sum_gff_length
#a subroutine that takes a gff file as input and returns the sum of all its features (no checking fo overlap or anything fancy)

sub sum_gff_length{
	use strict;
	use warnings;
	
	my ($file) = shift (@_) or die ("\n\nError in sum_gff_length: No file provided\n\n");
	my @fileData = file_to_array($file);
	my $length = 0;
	foreach my $line (@fileData){
		if ($line =~ /^\s*$/){next;}
		my @splitLine = split (/\t/, $line);
		unless (defined $splitLine[4]){
			print "UNDEFINED: $line"
		}
		$length = $length + $splitLine[4] - $splitLine[3] + 1;
	}
	
	return $length
}
#end sum_gff_length
###############################################################################

###############################################################################
#BEGIN find_orfs
#takes a nucleotide sequence, sequence name, strand (1 for fwd, 0 for rev) and minimum ORF length in (bp)
#returns a list of all orfs in gff format

sub find_orfs {
	use warnings;
	use strict;
	
	my ($sequence, $seqName, $minLength, $forwardStrand) = @_;
	my $gffString;
	
	unless ($forwardStrand){
		$sequence = revcom($sequence);
	}
	
	for (my $frame = 1; $frame <=3; ++$frame){
		my $startCodonPos;
		my $stopCodonPos;
		for (my $i = $frame - 1; $i < length($sequence) - 2; $i+=3){
			my $aa = codon2aa(substr($sequence, $i, 3));
			if ($aa eq 'M'){
				unless ($startCodonPos){
					$startCodonPos = $i + 1;
				}
			}elsif ($aa eq '*'){
				$stopCodonPos = $i+3;
				if ($startCodonPos){
					my $ORFlength = $stopCodonPos-$startCodonPos+1;
					if ($ORFlength >= $minLength){
						if ($forwardStrand){
							$gffString .= "$seqName\tfind_orfs\tORF\t$startCodonPos\t$stopCodonPos\t$ORFlength\t+\t$frame\t.\n"
						}else{
							my $revStart = length($sequence) - $stopCodonPos + 1;
							my $revStop = length($sequence) - $startCodonPos + 1;
							$gffString .= "$seqName\tfind_orfs\tORF\t$revStart\t$revStop\t$ORFlength\t-\t-$frame\t.\n"
						}
					}
				}
				$startCodonPos = undef;
				$stopCodonPos = undef;
			}
		}
		if ($startCodonPos){
			my $ORFlength = length($sequence)-$startCodonPos+1;
			
			if ($ORFlength >= $minLength){ #check for ORFs in progress (possibly interupted by break in circular chromosome
				if ($forwardStrand){
					$gffString .= "$seqName\tfind_orfs\tpartialORF\t$startCodonPos\t".length($sequence)."\t$ORFlength\t+\t$frame\t.\n"
				}else{
					my $revStart = 1;
					my $revStop = length($sequence) - $startCodonPos + 1;
					$gffString .= "$seqName\tfind_orfs\tpartialORF\t$revStart\t$revStop\t$ORFlength\t-\t-$frame\t.\n"
				}

			}
		}
	} 
	return $gffString;
}


###############################################################################
#BEGIN find_orfs_both_strands
#takes a sequence, seqname, and minimumLength and calls find_orfs for both strands
#returns concatenated output in gff format

sub find_orfs_both_strands{
	use warnings;
	use strict;
	
	my ($sequence, $seqName, $minLength) = @_;
	my $gffString;
	
	$gffString .= find_orfs($sequence, $seqName, $minLength, 1);
	$gffString .= find_orfs($sequence, $seqName, $minLength, 0);
	
	return $gffString;
}


###############################################################################
###############################################################################
#BEGIN genbank2fasta
#A subroutine that takes a genbank file and outputs the corresponding sequence in fasta format

sub genbank2fasta{
	use warnings;
	use strict;
	use Bio::SeqIO;
	
	my $gb_file = shift(@_) or die ("\n\nNo fasta file name provided to subroutine genbank2fasta\n\n");
	my $seqio_obj = Bio::SeqIO->new(-file => $gb_file);
	my $seq_obj = $seqio_obj->next_seq;
	
	my $header = ">".$seq_obj->id."|".$seq_obj->desc;
	my $sequence = $seq_obj->seq;
	
	return ($header, $sequence);
}
###############################################################################

###############################################################################
#BEGIN hash2fasta
#A subroutine that takes a hash of seq names (keys) and sequences (values) and exports it as a fasta file

sub hash2fasta{
	use warnings;
	use strict;

	my $usage = "\nERROR in hash2fasta: invalid arguments\n\n";
	my $hashRef = shift (@_) or die ($usage);
	my %hash = %$hashRef;
	my $file = shift (@_) or die ($usage);
	
	my $FH = open_output($file);
	
	foreach my $header (sort keys %hash){
		print $FH ">$header\n$hash{$header}\n";
	}
	
	close $FH;

}
###############################################################################

###############################################################################
#BEGIN align_fasta_with_muscle
#A subroutine that takes a fasta file and calls muscle for a nucleotide alignment. First argument is input. Optional second argument is output. Appends muscleAligned.fas to file name if output file not specified

sub align_fasta_with_muscle{
	use warnings;
	use strict;

	my $usage = "\nERROR in align_fasta_with_muscle: invalid arguments\n\n";
	my $inputFile = shift (@_) or die($usage);
	my $outputFile; 
	$outputFile = shift (@_) or $outputFile = "$inputFile\.muscleAligned.fas";
	
	system ("muscle -in $inputFile -out $outputFile -quiet");

}
###############################################################################
###############################################################################
#BEGIN align_fasta_with_mafft
#A subroutine that takes a fasta file and calls maaft for a nucleotide alignment. First argument is input. Optional second argument is output filename. Optional third argument is options string. Appends mafftAligned.fas to file name if output file not specified. If third optional argument is specified but not outputfilename, you must precede them with null arguments ("").

sub align_fasta_with_mafft{
	use warnings;
	use strict;
	

	my $usage = "\nERROR in align_fasta_with_muscle: invalid arguments\n\n";
	my $inputFile = shift (@_) or die($usage);
	my $outputFile; 
	$outputFile = shift (@_);
	$outputFile or $outputFile = "$inputFile\.maaftAligned.fas";
	my $options = shift;

	my $mafft;
	if ($options){
		$mafft = "mafft $options";
	}else{
		$mafft = "mafft";
	}
		
	system ("$mafft --quiet $inputFile > $outputFile");

}
###############################################################################
###############################################################################


#BEGIN align_fasta_with_muscle_codon
#A subroutine that takes a nucleotide fasta file, translates it, calls muscle for a protein alignment, then converts the aligned output back into nucleotide data. First argument is input. Optional second argument is output. Appends muscleAlignedCodon.fas to file name if output file not specified

sub align_fasta_with_muscle_codon{
	use warnings;
	use strict;

	my $usage = "\nERROR in align_fasta_with_muscle_codon: invalid arguments\n\n";
	my $inputFile = shift (@_) or die($usage);
	my $outputFile; 
	$outputFile = shift (@_) or $outputFile = "$inputFile\.muscleAlignedCodon.fas";
	
	my %inputFastaHash = arrays2hash(get_fasta_names_and_seqs($inputFile));
	
	my $FH = open_output ("$inputFile\.tempTranslationForMuscle.fas");
	
	foreach my $header (sort keys %inputFastaHash){
		my $aaSeq = dna2peptide($inputFastaHash{$header});
		if (substr($aaSeq, -1, 1) eq '*'){
			$aaSeq = substr($aaSeq, 0, -1);
			print STDERR "\nWARNING: Removing stop codon from 3' end of $header\n";
		}
		if ($aaSeq =~ /\*/){
			die ("\nERROR in align_fasta_with_muscle_codon: Sequence cannot contain internal stop codons\n\n");
		}
		$aaSeq =~ s/\?/X/g; 
		print $FH ">$header\n$aaSeq\n";
	}
	close $FH;
		
	system ("muscle -in $inputFile\.tempTranslationForMuscle.fas -out $inputFile\.tempTranslationForMuscle.aligned.fas -quiet");
	
	my %alignedHash = arrays2hash(get_fasta_names_and_seqs("$inputFile\.tempTranslationForMuscle.aligned.fas"));
	
	system ("rm -f $inputFile\.tempTranslationForMuscle*");
	
	my $FH_FINAL = open_output ($outputFile);
	foreach my $header (sort keys %inputFastaHash){
		my $seq;
		my $nucPos = 0;
		for (my $i = 0; $i < length ($alignedHash{$header}); ++$i){
			if (substr($alignedHash{$header}, $i, 1) eq '-'){
				$seq .= '---';
			}else{
				$seq .= substr($inputFastaHash{$header},$nucPos, 3);
				$nucPos += 3;
			}
		}
		print $FH_FINAL ">$header\n$seq\n"; 
	}
	close $FH_FINAL;	
}
###############################################################################




###############################################################################
###############################################################################
#BEGIN align_fasta_with_mafft_codon
#A subroutine that takes a nucleotide fasta file, translates it, calls maaft for a protein alignment, then converts the aligned output back into nucleotide data. First argument is input. Second argument is output. Third argument is optional string of maaft command line options

sub align_fasta_with_mafft_codon{
	use warnings;
	use strict;

	my $usage = "\nERROR in align_fasta_with_mafft_codon: invalid arguments\n\n";
	my $inputFile = shift (@_) or die($usage);
	my $outputFile = shift or die ($usage);
	my $options = shift;
	
	my %inputFastaHash = arrays2hash(get_fasta_names_and_seqs($inputFile));
	
	my $FH = open_output ("$inputFile\.tempTranslationForMAFFT.fas");
	
	foreach my $header (sort keys %inputFastaHash){
		my $aaSeq = dna2peptide($inputFastaHash{$header});
		if (substr($aaSeq, -1, 1) eq '*'){
			$aaSeq = substr($aaSeq, 0, -1);
			print STDERR "\nWARNING: Removing stop codon from 3' end of $header\n";
		}
		if ($aaSeq =~ /\*/){
			die ("\nERROR: in align_fasta_with_mafft_codon: Sequence contains internal stop codons\n\n");
		}
		$aaSeq =~ s/\?/X/g; 
		print $FH ">$header\n$aaSeq\n";
	}
	close $FH;
	
	my $mafft;
	if ($options){
		$mafft = "mafft $options";
	}else{
		$mafft = "mafft";
	}
		
	#system ("$mafft --quiet $inputFile\.tempTranslationForMAFFT.fas > $inputFile\.tempTranslationForMAFFT.aligned.fas");
	system ("einsi --quiet $inputFile\.tempTranslationForMAFFT.fas > $inputFile\.tempTranslationForMAFFT.aligned.fas");
	
	my %alignedHash = arrays2hash(get_fasta_names_and_seqs("$inputFile\.tempTranslationForMAFFT.aligned.fas"));
	
	system ("rm -f $inputFile\.tempTranslationForMAFFT*");
	
	my $FH_FINAL = open_output ($outputFile);
	foreach my $header (sort keys %inputFastaHash){
		my $seq;
		my $nucPos = 0;
		for (my $i = 0; $i < length ($alignedHash{$header}); ++$i){
			if (substr($alignedHash{$header}, $i, 1) eq '-'){
				$seq .= '---';
			}else{
				$seq .= substr($inputFastaHash{$header},$nucPos, 3);
				$nucPos += 3;
			}
		}
		print $FH_FINAL ">$header\n$seq\n"; 
	}
	close $FH_FINAL;	
}
###############################################################################


###############################################################################
#BEGIN download_genbank_file
#A subroutine that takes a GenBank accession number and a directory and downloads the corresponding genbank file to the directory 

sub download_genbank_file{
	use warnings;
	use strict;
	use Bio::DB::GenBank;
	use Bio::SeqIO;

	my $usage = "\nERROR in download_genbank_file: invalid arguments\n\n";
	my $accession = shift (@_) or die($usage);
	my $directory = shift (@_) or die($usage);
	
	unless (substr($directory, -1, 1) eq '/'){
		$directory .= '/';
	} 

	my $db_obj = Bio::DB::GenBank->new;
	my $seq_obj = $db_obj->get_Seq_by_acc($accession);
	my $seqio_obj = Bio::SeqIO->new(-file => ">$directory$accession\.gb", -format => 'genbank' );
	
	$seqio_obj->write_seq($seq_obj);

}

###############################################################################
###############################################################################
#BEGIN download_genbank_file_multiple
#A subroutine that takes a ref to an array of GenBank accession numbers and a directory and downloads the corresponding genbank files to the directory 

sub download_genbank_file_multiple{
	use warnings;
	use strict;
	use Bio::DB::GenBank;
	use Bio::SeqIO;

	my $usage = "\nERROR in download_genbank_file_multiple: invalid arguments\n\n";
	my $accessionArrayRef = shift (@_) or die($usage);
	my @accessions = @$accessionArrayRef;
	my $directory = shift (@_) or die($usage);
	
	unless (substr($directory, -1, 1) eq '/'){
		$directory .= '/';
	} 

	my $db_obj = Bio::DB::GenBank->new;
	
	#note that the square brackets around @accessiona are necessary to reference array
	my $stream_obj = $db_obj->get_Stream_by_acc([@accessions]);
	
	my $count = 0;
	while (my $seq_obj = $stream_obj->next_seq()){
		if ($seq_obj->length()>0){ #only proceed if a sequence was successfully downloaded
			my $seqio_obj = Bio::SeqIO->new(-file => ">$directory$accessions[$count]\.gb", -format => "genbank" );
			$seqio_obj->write_seq($seq_obj);			
		} else {
			print STDERR "WARNING. FAILED DOWNLOAD. $accessions[$count]";
			++$count;
		}
		++$count;
	}
}

###############################################################################
###############################################################################
#BEGIN fasta2phy
#A subroutine that takes a fasta file and converts it to a phylip file. Second (optional) argument is a delimeter between name and sequence.

sub fasta2phy {
	my $fasta = shift @_ or die ("\nERROR: No fasta file name provided to fasta2phy\n\n");
	my $delim;
	$delim = shift @_ or $delim = " ";
	my %fastaHash = arrays2hash(get_fasta_names_and_seqs($fasta));
	
	my @seqNames = sort keys %fastaHash;
	my $numSeqs = scalar @seqNames;
	my $length = length ($fastaHash{$seqNames[0]});
	
	foreach my $name (@seqNames){
		length($fastaHash{$name}) != $length && die ("\nERROR in fasta2phy: Not all sequences are of equal length\n\n");
		$name =~ /\s+/ && die ("\nERROR in fasta2phy: sequence names cannot contain whitespace\n\n");
	}
	
	my $phy = "$numSeqs $length\n";
	foreach my $name (@seqNames){
		$phy .= "$name$delim$fastaHash{$name}\n";
	}
	
	return $phy;
	
}
###############################################################################
###############################################################################
#BEGIN phy2fasta
#A subroutine that takes a phylip alignment (sequential) and converts it to a fasta file

sub phy2fasta {
	my $phylip = shift @_ or die ("\nERROR: No phylip file name provided to phy2fasta\n\n");
	
	my @phyLines = file_to_array($phylip);
	my $fasta;
	for (my $i = 1; $i < scalar (@phyLines); ++$i){
		chomp $phyLines[$i];
		my ($name, $seq) = split (/\s+/, $phyLines[$i]);
		$fasta .= ">$name\n$seq\n";
	}
		
	return $fasta;
	
}

###############################################################################
#BEGIN fasta2hash
#A subroutine that takes a fasta file and returns a hash with headers as keys and sequences as values

sub fasta2hash {
	my $fasta = shift @_ or die ("\nERROR: No fasta file name provided to fasta2hash\n\n");
	my %fastaHash = arrays2hash (get_fasta_names_and_seqs($fasta));
	return %fastaHash;
}


1;