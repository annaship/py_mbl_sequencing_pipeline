
import sys, os, stat



class IlluminaFiltering:
    """reads an Illumina fastq file and outputs a new filtered file
         Filtering can be on chastity, Ns, mate exists, B-tails, etc.
         Based partially on minoche_filtering_pipeline.pl, but with more controls
         /bioware/illumina_scripts/minoche_filtering_pipeline.pl
         
    Usage:  illumina_filtering.pl -chastity -N -Btail -len min_length -trim trim_length -clip clip_bases -unmated mate.fastq -in input.fastq -failed failed.fastq -out output.fastq -outmate outmate.fastq

      ex:  illumina_filtering.pl -chastity -in PhiX_TTAGGC_L008_R1_001.fastq -out PhiX_TTAGGC_L008_R1_001.chaste.fastq 
           illumina_filtering.pl -chastity -N -in PhiX_TTAGGC_L008_R1_001.fastq -failed PhiX_TTAGGC_L008_R1_001.failed.fastq -out PhiX_TTAGGC_L008_R1_001.chaste.fastq 
           illumina_filtering.pl -unmated PhiX_TTAGGC_L008_R2_001.chaste.N.fastq -outmate PhiX_TTAGGC_L008_R2_001.chaste.N.mated.fastq -in PhiX_TTAGGC_L008_R1_001.chaste.N.fastq -out PhiX_TTAGGC_L008_R1_001.chaste.N.mated.fastq 

 Options:  
           -in          input fastq file
           -out         output fastq file of reads passing quality filters
           -failed      output fastq file of reads failing quality filters
           -unmated     input mate fastq file to compare for removing reads that no longer have a mate
           -outmate     output fastq file for mated reads from -unmated input file

           -chastity    remove reads that do not pass the Illumina chastity filter
           -N           remove reads that contain an ambiguous base
           -Btail       remove reads that contain a Btail
           -unmated     remove reads that no longer have a mate in the mate fastq file
           -len         remove reads that are shorter than minimum length
           -f50         remove reads that have >=34 qual scores < 30 in the first 50 nt (>2/3s)
           -trim        trim remaining reads to trim length
           -clip        clip N bases from the start of remaining reads (for R2)
           -Nx          ignore ambiguous base at position X (seems to be in the v6 primer)
    """
    Name = "IlluminaFiltering"
    def __init__(self, run_object, infile=None, outfile=None, chastity=False, filter_first50=False, filter_Ns=False, filter_Nx=0, trim=0, clip=0, btail=False, length=0, failed=None ):
        self.infile         = infile
        self.outfile        = outfile
        self.failed         = failed
        self.filter_chaste  = chastity
        self.filter_first50 = filter_first50
        self.filter_Ns      = filter_Ns
        self.filter_Nx      = filter_Nx
        self.trim_length    = trim
        self.clip_length    = clip
        self.filter_Btails  = btail
        self.filter_length  = length
        self.runobj         = run_object
        self.indir          = self.runobj.input_dir
        self.outdir         = os.path.join(self.runobj.output_dir,"trim")
    
    def chastity_filter(self, file):
        from Bio import SeqIO
        
        
        #illumina_filtering -chastity -N -f50 -in $1.fastq -out $1.chaste.N.f50.fastq -failed $1.chaste.N.f50.failed.fastq
        #   -chastity    remove reads that do not pass the Illumina chastity filter
        #   -N           remove reads that contain an ambiguous base
        #   -f50         remove reads that have >=34 qual scores < 30 in the first 50 nt (>2/3s)
        # # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
        #       0           1               2           3      4     5         6  -  6        7             8                   9
        infile = os.path.join(self.indir,file)
        filebase        = file.split('/')[1].split('.')[0]
        output_handle   = open(os.path.join(self.outdir, filebase+".chaste.N.f50.fastq"), "w")
        fail_handle     = open(os.path.join(self.outdir, filebase+".chaste.N.f50.failed"), "w")
        in_handle = open(infile, "rU")
        count_of_unchaste = 0
        for record in SeqIO.parse(in_handle, "fastq-sanger") :
            #print record.id
            print record.description
            seq = record.seq
            desc_items = record.description.split(':')
            
            if desc_items[7] == 'Y':
                count_of_unchaste += 1
                count = SeqIO.write(record, fail_handle, "fastq-sanger")
                continue
            
            
            if self.filter_Ns:
                print '',seq[self.filter_Nx-1:self.filter_Nx]
                countN = seq.count('N')
                if countN > 1 or (countN == 1 and seq[self.filter_Nx-1:self.filter_Nx] != 'N'):
                    count = SeqIO.write(record, fail_handle, "fastq-sanger")
                    continue
                    
            if self.filter_length:
                if len(seq) < self.filter_length:
                    count = SeqIO.write(record, fail_handle, "fastq-sanger")
                    continue
                    
            if self.filter_first50:
                record.quality
                
                
            if self.clip_length:
                pass
                
            if self.trim_length:
                pass
                
            # write out the clean reads
            SeqIO.write(record, output_handle, "fastq-sanger")
                
        in_handle.close()
        
        
    def btails_filter(self, infile):
        """
        """
        #illumina_filtering -Btail -in $1.chaste.N.f50.fastq -out $1.chaste.N.f50.Btails.fastq -failed $1.chaste.N.f50.Btails.failed
        #   -Btail       remove reads that contain a Btail
        in_filename = infile
        
        print "Filtering Btails"
        report_filename = self.outfile + ".btail.report";
        error_filename = self.outfile + ".btail.errors";
        print "Running Galaxy fastq_trimmer_by_quality to remove B tails...this can take several hours on a large file...\n";
        print "Output from fastq_trimmer_by_quality:\n >> $report_filename";
        btail_cmd = "/bioware/galaxy/tools/fastq/fastq_trimmer_by_quality.py -s 1 -t 1 -e 53 -a mean -x 0 -c '>=' -q 3 "+in_filename+" "+out_filename+" >> "+report_filename+" 2>"+error_filename
        
        subprocess.call(btail_cmd,shell=True)
    
    
        
        
    def length_filter(self):
        """
        """
        #illumina_filtering -len 75 -in $1.chaste.N.f50.Btails.fastq -out $1.chaste.N.f50.Btails.len75.fastq -failed $1.chaste.N.f50.Btails.len75.failed.fastq
        #   -len         remove reads that are shorter than minimum length
        #   -clip        clip N bases from the start of remaining reads (for R2)
        
#######################################
#
# Set up usage statement
#
#######################################
# my $script_help = "
#  illumina_filtering.pl - reads an Illumina fastq file and outputs a new filtered file
#                          Filtering can be on chastity, Ns, mate exists, B-tails, etc.
#                          Based partially on minoche_filtering_pipeline.pl, but with more controls
#                          /bioware/illumina_scripts/minoche_filtering_pipeline.pl
# \n";
# 
# my $usage = "
#    Usage:  illumina_filtering.pl -chastity -N -Btail -len min_length -trim trim_length -clip clip_bases -unmated mate.fastq -in input.fastq -failed failed.fastq -out output.fastq -outmate outmate.fastq
# 
#       ex:  illumina_filtering.pl -chastity -in PhiX_TTAGGC_L008_R1_001.fastq -out PhiX_TTAGGC_L008_R1_001.chaste.fastq 
#            illumina_filtering.pl -chastity -N -in PhiX_TTAGGC_L008_R1_001.fastq -failed PhiX_TTAGGC_L008_R1_001.failed.fastq -out PhiX_TTAGGC_L008_R1_001.chaste.fastq 
#            illumina_filtering.pl -unmated PhiX_TTAGGC_L008_R2_001.chaste.N.fastq -outmate PhiX_TTAGGC_L008_R2_001.chaste.N.mated.fastq -in PhiX_TTAGGC_L008_R1_001.chaste.N.fastq -out PhiX_TTAGGC_L008_R1_001.chaste.N.mated.fastq 
# 
#  Options:  
#            -in          input fastq file
#            -out         output fastq file of reads passing quality filters
#            -failed      output fastq file of reads failing quality filters
#            -unmated     input mate fastq file to compare for removing reads that no longer have a mate
#            -outmate     output fastq file for mated reads from -unmated input file
# 
#            -chastity    remove reads that do not pass the Illumina chastity filter
#            -N           remove reads that contain an ambiguous base
#            -Btail       remove reads that contain a Btail
#            -unmated     remove reads that no longer have a mate in the mate fastq file
#            -len         remove reads that are shorter than minimum length
#            -f50         remove reads that have >=34 qual scores < 30 in the first 50 nt (>2/3s)
#            -trim        trim remaining reads to trim length
#            -clip        clip N bases from the start of remaining reads (for R2)
#            -Nx          ignore ambiguous base at position X (seems to be in the v6 primer)
# \n";
# 
# #######################################
# #
# # Definition statements
# #
# #######################################
# #Commandline parsing
# my $arg_count = 0;
# my $verbose = 0;
# my $self_cmd = join(" ", $0, @ARGV);
# 
# #Runtime variables
# my $in_filename;
# my $out_filename;
# my $mate_filename;
# my $outmate_filename;
# my $gzipped_input = 0;
# my $filter_chaste = 0;
# my $filter_Ns = 0;
# my $filter_Nx = 0;
# my $filter_Btails = 0;
# my $filter_unmated = 0;
# my $failed_fastq = 0;
# my $filter_length = -1;
# my $filter_first50 = 0;
# my $first50 = 50;
# my $first50_maxQ = 30;
# my $first50_maxQ_count = 34;
# my $trim_length = 0;
# my $clip_length = 0;
# 
# my $count_of_unchaste = 0;
# my $count_of_Ns = 0;
# my $count_of_shorts = 0;
# my $count_of_unmated = 0;
# my $count_of_first50 = 0;
# my $count_of_trimmed = 0;
# 
# my $log_filename = $0;
# $log_filename =~ s/^.*\///;
# $log_filename = "./" . $log_filename . ".log";
# 
# #######################################
# #
# # Test for commandline arguments
# #
# #######################################
# 
# if (! $ARGV[0] ) 
# {
# 	print $script_help;
# 	print $usage;
# 	exit -1;
# } 
# 
# 
# while ((scalar @ARGV > 0) && ($ARGV[0] =~ /^-/)) 
# {
# 	if ($ARGV[0] =~ /-h/) 
# 	{
# 		print $script_help;
# 		print $usage;
# 		exit 0;
# 	} elsif ($ARGV[0] eq "-in") {
# 		shift @ARGV;
# 		$in_filename = shift @ARGV;
# 	} elsif ($ARGV[0] eq "-out") {
# 		shift @ARGV;
# 		$out_filename = shift @ARGV;
# 	} elsif ($ARGV[0] eq "-outmate") {
# 		shift @ARGV;
# 		$outmate_filename = shift @ARGV;
# 	} elsif ($ARGV[0] eq "-chastity") {
# 		shift @ARGV;
# 		$filter_chaste = 1;
# 	} elsif ($ARGV[0] eq "-f50") {
# 		shift @ARGV;
# 		$filter_first50 = 1;
# 	} elsif ($ARGV[0] eq "-N") {
# 		shift @ARGV;
# 		$filter_Ns = 1;
# 	} elsif ($ARGV[0] eq "-Nx") {
# 		shift @ARGV;
# 		$filter_Nx = shift @ARGV;
# 	} elsif ($ARGV[0] eq "-Btail") {
# 		shift @ARGV;
# 		$filter_Btails = 1;
# 	} elsif ($ARGV[0] eq "-len") {
# 		shift @ARGV;
# 		$filter_length = shift @ARGV;
#         if ($filter_length !~ /^[0-9]+$/)
#         {
#             print "-len argument must be a positive integer value. Exiting.\n"; exit;
#         }
# 	} elsif ($ARGV[0] eq "-trim") {
# 		shift @ARGV;
# 		$trim_length = shift @ARGV;
#         if ($trim_length !~ /^[0-9]+$/)
#         {
#             print "-trim argument must be a positive integer value. Exiting.\n"; exit;
#         }
#         #print "Trim to length: $trim_length\n";
# 	} elsif ($ARGV[0] eq "-clip") {
# 		shift @ARGV;
# 		$clip_length = shift @ARGV;
#         if ($clip_length !~ /^[0-9]+$/)
#         {
#             print "-clip argument must be a positive integer value. Exiting.\n"; exit;
#         }
# 	} elsif ($ARGV[0] eq "-unmated") {
# 		shift @ARGV;
# 		$filter_unmated = 1;
#         $mate_filename = shift @ARGV;
# 	} elsif ($ARGV[0] eq "-failed") {
# 		shift @ARGV;
# 		$failed_fastq = shift @ARGV;
# 	} elsif ($ARGV[0] eq "-v") {
# 		$verbose = 1;
# 		shift @ARGV;
# 	} elsif ($ARGV[0] =~ /^-/) { #unknown parameter, just get rid of it
# 		print "Unknown commandline flag \"$ARGV[0]\".\n";
# 		print $usage;
# 		exit -1;
# 	}
# }
# 
# 
# #######################################
# #
# # Parse commandline arguments, ARGV
# #
# #######################################
# 
# if (scalar @ARGV != $arg_count) 
# {
# 	print "Incorrect number of arguments.\n";
# 	print "$usage\n";
# 	exit;
# } 
# 
# if ( (! $in_filename) || (! $out_filename) )
# {
# 	print "Must specify both input and output fastq filenames.\n";
# 	print "$usage\n";
# 	exit;
# } 
# 
# if ( ($filter_Btails) && (($filter_chaste) || ($filter_Ns) || ($filter_unmated)) )
# {
# 	print "Filtering Btails should be done independently of all other filtering options, sorry.\n";
# 	exit -1;
# }
#     
# if ( ($filter_unmated) && (($filter_chaste) || ($filter_Ns) || ($filter_Btails)) )
# {
# 	print "Filtering unmated should be done independently of all other filtering options, sorry.\n";
# 	exit -1;
# }
#     
# 
# #Test validity of commandline arguments
# if ( ($in_filename ne "stdin") && (! -f $in_filename) ) 
# {
# 	print "Unable to locate input fastq file: $in_filename.\n";
# 	exit -1;
# }
# 
# if ( ($mate_filename) && (! -f $mate_filename) ) 
# {
# 	print "Unable to locate input mate fastq file: $mate_filename.\n";
# 	exit -1;
# }
# 
# open(LOG, ">>$log_filename")  || warn "Unable to open log file, $log_filename, for writing.  Exiting...\n";
# print LOG "$self_cmd\n";
# if ( ($in_filename =~ /\.gz$/) || ($in_filename =~ /\.gzip$/) )
# {
#     $gzipped_input = 1;
#     my $gunzip_error = system("gunzip $in_filename");
#     if ($gunzip_error) {print "Unable to gunzip $in_filename.  Exiting\n\n"; exit -1;}
#     if ($in_filename =~ /\.gzip$/) 
#     {
#         $in_filename =~ s/\.gzip$//;
#     } else {
#         $in_filename =~ s/\.gz$//;
#     }
# }
# 
# #######################################
# #
# # Remove Btails using fastq_trimmer_by_quality.py from galaxy 
# #
# #######################################
#     #
#     # Usage: fastq_trimmer_by_quality.py [options] input_file output_file
#     #
#     # Options:
#     #   -h, --help            show this help message and exit
#     #   -f FORMAT, --format=FORMAT
#     #                          FASTQ variant type
#     #   -s WINDOW_SIZE, --window_size=WINDOW_SIZE
#     #                          Window size
#     #   -t WINDOW_STEP, --window_step=WINDOW_STEP
#     #                          Window step
#     #   -e TRIM_ENDS, --trim_ends=TRIM_ENDS
#     #                          Ends to Trim
#     #   -a AGGREGATION_ACTION, --aggregation_action=AGGREGATION_ACTION
#     #                          Aggregate action for window
#     #   -x EXCLUDE_COUNT, --exclude_count=EXCLUDE_COUNT
#     #                          Maximum number of bases to exclude from the window
#     #                          during aggregation
#     #   -c SCORE_COMPARISON, --score_comparison=SCORE_COMPARISON
#     #                          Keep read when aggregate score is
#     #   -q QUALITY_SCORE, --quality_score=QUALITY_SCORE
#     #                          Quality Score
#     #   -k, --keep_zero_length
#     #                          Keep reads with zero length
#     #
# 
# if ($filter_Btails)
# {
#     print "Filtering Btails\n";
#     my $report_filename = $out_filename . "btail.report";
#     my $error_filename = $out_filename . "btail.errors";
#     open my $rfh, ">>$report_filename";
#     print $rfh "Output from fastq_trimmer_by_quality:\n";
# 
#     print "Running Galaxy fastq_trimmer_by_quality to remove B tails...this can take several hours on a large file...\n";
#     system("python /bioware/galaxy/tools/fastq/fastq_trimmer_by_quality.py -s 1 -t 1 -e 53 -a mean -x 0 -c '>=' -q 3 $in_filename $out_filename >> $report_filename 2>$error_filename");
#     print $rfh "\n";
#     exit;
# }
# 
# 
# #######################################
# #
# # Perform all other filtering
# #
# #######################################
# 
# if ( ($filter_chaste) || ($filter_Ns) || ($filter_length > -1) || ($filter_first50) || ($trim_length) )
# {
#     print "Filtering reads\n";
#     #
#     # Open the files
#     #
#     
#     # Open the input fastq file
#     my $in;
#     if ($in_filename eq "stdin")
#     {   
#         $in = Bio::SeqIO->new( '-fh'=> \*STDIN, '-format'=> "fastq") || die("Could not read input fastq file from STDIN.  Exiting...\n");
#     } else {
#         $in = Bio::SeqIO->new( '-file'=> "<$in_filename", '-format'=> "fastq", -variant => 'sanger') || die("Could not read input fastq file: $in_filename.  Exiting...\n");
#     }   
#     
#     # Open the output fastq
#     my $out = Bio::SeqIO->new( '-file'=> ">$out_filename", '-format'=> "fastq", -variant => 'sanger') || die("Unable to write to output file: $out_filename.  Exiting...\n");
#     
#     # Open the failed reads output fastq
#     my $fail;
#     if ($failed_fastq)
#     {
#         $fail = Bio::SeqIO->new( '-file'=> ">$failed_fastq", '-format'=> "fastq", -variant => 'sanger') || die("Unable to write to output file: $failed_fastq.  Exiting...\n");
#     }
#     
#     #
#     # Step through the sequences and print only chaste reads
#     #
#     while (my $seq_object = $in->next_seq) 
#     {
#        # load the data
#        #my $id = $seq_object->id();
#     
#        # Filter chaste reads
#        if ($filter_chaste) 
#        {
#            my $desc = $seq_object->desc();
#            if ($desc =~ /:Y:/) 
#            {
#                if ($failed_fastq) { $fail->write_seq($seq_object); }
#                $count_of_unchaste++;
#                next;
#            }
#        }
#     
#        my $seq = $seq_object->seq();
# 
#        # Filter reads with ambiguous bases
#        if ($filter_Ns) 
#         {
#             my $tmpseq = $seq;
#             my $countN = $tmpseq =~ tr/N/N/;
#             # is there only one N at the designated ignore position
#             if ( ($countN > 1) || ( ( $countN == 1) && (substr($seq, $filter_Nx - 1, 1) ne 'N') ) )
#             {        
#                 if ($failed_fastq) { $fail->write_seq($seq_object); }
#                 $count_of_Ns++;
#                 next;
#             }
#             $tmpseq='';
#         }
# 
#         # Filter reads below minimum length
#         if ($filter_length) 
#         {
#             if (length($seq) < $filter_length)
#             {
#                 if ($failed_fastq) { $fail->write_seq($seq_object); }
#                 $count_of_shorts++;
#                 next;
#             }
#         }
#     
#         # Filter reads below first 50 base quality
#         if ($filter_first50) 
#         {
#             my $quals = $seq_object->subqual(1, $first50);
#             my $count_lt30 = 0;
# 
#             # Step through the qual values and check if they are less than 30 (maxQ)
#             foreach my $q (@$quals)
#             {
#                 if ($q < $first50_maxQ) {$count_lt30++;}
#             }
# 
#             if ($count_lt30 >= $first50_maxQ_count)
#             {
#                 if ($failed_fastq) { $fail->write_seq($seq_object); }
#                 $count_of_first50++;
#                 next;
#             }
#                 #print "Pass $count_lt30: " . join(" ", @$quals) . "\n";
#         }
# 
#         # Trim initial bases -- remove first 10 bases from read 2
#         if ($clip_length)
#         {
#             $seq = substr($seq, $clip_length, length($seq) );
#             $seq_object->seq($seq);
# 
#             $seq_object->qual( $seq_object->subqual($clip_length+1, length($seq)) );
#             $count_of_trimmed++;
#         }
#     
#         # Trim to max length -- read 2 trim to 90.
#         if ($trim_length)
#         {
#             if (length($seq) > $trim_length)
#             {
#                 my $seq = substr($seq, 0, $trim_length);
#                 $seq_object->seq($seq);
# 
#                 $seq_object->qual( $seq_object->subqual(1, $trim_length) );
#                 $count_of_trimmed++;
#             }
#         }
#     
#         # Write out the clean reads
#     	$out->write_seq($seq_object);
#     }
# }
# 
# #######################################
# #
# # Filter to mated only
# #
# #######################################
# if ($filter_unmated)
# {
#     print "Filtering unmated\n";
#     
#     # Open the input fastq files
#     my $in = Bio::SeqIO->new( '-file'=> "<$in_filename", '-format'=> "fastq", -variant => 'sanger') || die("Could not read input fastq file: $in_filename.  Exiting...\n");
#     my $mate = Bio::SeqIO->new( '-file'=> "<$mate_filename", '-format'=> "fastq", -variant => 'sanger') || die("Could not read input fastq file: $mate_filename.  Exiting...\n");
# 
#     # Open the output fastq files
#     my $out = Bio::SeqIO->new( '-file'=> ">$out_filename", '-format'=> "fastq", -variant => 'sanger') || die("Unable to write to output file: $out_filename.  Exiting...\n");
#     my $outmate = Bio::SeqIO->new( '-file'=> ">$outmate_filename", '-format'=> "fastq", -variant => 'sanger') || die("Unable to write to output file: $outmate_filename.  Exiting...\n");
# 
#     # Load up the mate
#     my %still_in_mate;
#     while (my $seq_object = $mate->next_seq) 
#     {
#         # load the data
#     	my $id = $seq_object->id();
#         if ($id =~ /2308:15139:199300/) {print "Mate is: $id\n";}
#         $still_in_mate{$id} = $seq_object;
#     }
# 
#     # Step through infile and export kepts
#     while (my $seq_object = $in->next_seq) 
#     {
#         # Check each ID
#     	my $id = $seq_object->id();
# 
#         if ($id =~ /2308:15139:199300/) {print "Found in source\n";}
#         # Print the ones that are still in the mate file
#         if (exists $still_in_mate{$id}) 
#         { 
#             if ($id =~ /2308:15139:199300/) {print "Found in hash\n";}
#             $out->write_seq($seq_object); 
#             $outmate->write_seq($still_in_mate{$id}); 
#             delete $still_in_mate{$id};  #free up some memory
#         }
#     }
# }
#     
# if ($gzipped_input) 
# {
#     my $gzip_error = system("gzip $in_filename");
#     if ($gzip_error) {print "Unable to gzip $in_filename. File has not been recompressed, sorry. \n\n"; exit -1;}
# }
