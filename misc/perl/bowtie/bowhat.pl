#!/usr/bin/env perl
#
# Kyle McChesney
# bowhat.pl
# 
# Do tophat on bunch of files
# just the alignment
use warnings;
use strict;
use File::Slurp;
use Data::Printer;
use Getopt::Long;
use feature 'say';
use FileDo;

# config options hash
my %opts;
my $opts_ref;

GetOptions(
    "reads-directory=s"   => \$opts{'reads_directory'},
    "map-directory=s"     => \$opts{'map_directory'},
    "map-basename=s"      => \$opts{'map_basename'},
    "results-directory=s" => \$opts{'results_directory'},
    "cores=i"       	  => \$opts{'cores'},
    "matched=i"           => \$opts{'matched'},
    "size-file=s"         => \$opts{'size_file'},
    "config=s"            => \$opts{'config_file'},
    "bowtie=i"			  => \$opts{'bowtie'},
    "tophat=i"            => \$opts{'tophat'},
    "gtf-file=s"          => \$opts{'gtf_file'}
) or die("malformed command line args \n");

# set params from config
if (defined($opts{'config_file'})) {
	$opts_ref = generate_config(\%opts);
}

my $ts = time;
$$opts_ref{'ts'} = $ts;
say "Begining run with the following options:";
p($opts_ref);

say "creating the results directory if needed";
unless(-d $$opts_ref{'results_directory'}) {
	say "mkdir ".$$opts_ref{'results_directory'};
}

say "creating results subdirectory";
say "mkdir".$$opts_ref{'results_directory'}."bowhat-out-".$$opts_ref{'ts'}."/";

if($$opts_ref{'matched'} == 1) {
	say "generating matched pairs";
	my @read_files_m = grep (/.*\.fastq/, read_dir($$opts_ref{'reads_directory'}));
	my %groups;
	foreach my $file (@read_files_m) {
		$file =~ m/([^_]*)(_\d)?\.fastq$/;
		my $basename = $1;
		
		if($2) {
			my $mate_num = $2;
			$mate_num =~ s/_//g;
			$groups{$basename}{$mate_num} = $$opts_ref{'reads_directory'}.$file; 
		} else {
			$groups{$basename}{0} = $$opts_ref{'reads_directory'}.$file;
		}
	}
	p(%groups);

	if($$opts_ref{'bowtie'} == 1) {
		say "preparing bowtie2 run";
		my $base_bowtie_command = "bowtie2 -p (?) -t --no-unal -x (?)(?) -s (?)";
		foreach my $group (keys(%groups)) {
			if(scalar(keys($groups{$group})) == 3) {
				my $full_command = $base_bowtie_command." -1 (?) -2 (?) -u (?)";
				my $binder = Binder->new(base_string => $full_command);
				$binder->bind(
					$$opts_ref{'cores'},
					$$opts_ref{'map_directory'},
					$$opts_ref{'map_basename'},
					$$opts_ref{'results_directory'}.$group.".bowtie.sam",
					$groups{$group}{1},
					$groups{$group}{2},
					$groups{$group}{0},
				);
				say "runing: ".$binder->bound_string();
				my $commnad = $binder->bound_string();
				say `$commnad`;
			}
		}
	}

	if($$opts_ref{'tophat'} == 1) {
		say "preparing tophat run";
		my $base_bowtie_command = "tophat2 -G (?) -p (?) -o (?) ";
		foreach my $group (keys(%groups)) {
			if(scalar(keys($groups{$group})) == 3) {
				my $full_command = $base_bowtie_command." -1 (?) -2 (?) -u (?)";
				my $binder = Binder->new(base_string => $full_command);
				$binder->bind(
					$$opts_ref{'cores'},
					$$opts_ref{'map_directory'},
					$$opts_ref{'map_basename'},
					$$opts_ref{'results_directory'}.$group.".bowtie.sam",
					$groups{$group}{1},
					$groups{$group}{2},
					$groups{$group}{0},
				);
				say "runing: ".$binder->bound_string();
				my $commnad = $binder->bound_string();
				say `$commnad`;
			}
		}

		#tophat2  -G $$opts{'gtf-file'} -p 5 -o $output_dir $$opts{'index-base'} $$opts{'reads-directory'}/$file_one $$opts{'reads-directory'}/$file_two
	}

}


sub generate_config {
	my $opts_hash = shift;
	my $config = $$opts_hash{'config_file'};
	my @lines = read_file($$opts_hash{'config_file'});

	foreach my $line (@lines) {
		
		# validate
		next if($line =~ m/^#/);
		my @arrayed = split(/\t/, $line);
		if(scalar(@arrayed) > 2) {
			say "@arrayed is malformed, skipping";
			next;
		}
		chomp($arrayed[1]);
		$$opts_hash{$arrayed[0]} = $arrayed[1];
	}
	return $opts_hash;
}
