#!/usr/bin/perl

#use strict;
#use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::SeqIO;
use Bio::EnsEMBL::TranscriptMapper;	

my $reg = "Bio::EnsEMBL::Registry";

print "Connect to DB\n";
$reg->load_registry_from_db(
    -host , "ensembldb.ensembl.org",
    -user , "anonymous"
);

my $comparadb= new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host   => "ensembldb.ensembl.org",
    -port   => 5306,
    -user   => "anonymous",
    -dbname => 'ensembl_compara_51',
);

print "Get Adaptors\n"; #(Species, Type of DB, Type of Object)
my $hg_adaptor = $reg->get_adaptor("Human","Core","Gene"); 
my $familyA = $reg->get_adaptor("Multi", "Compara", "Family");
my $seqmemberA = $reg->get_adaptor("Multi", "Compara", "SeqMember");
my $slice_adaptor = $reg->get_adaptor("Human", "Core", "Slice");
my $transcript_adaptor = $reg->get_adaptor("Human", "Core", "Transcript");
my $vf_adaptor = $reg->get_adaptor("Human", "Variation", "Variationfeature");
my $pa = $reg->get_adaptor("Human","Variation","Phenotype");
my $gene_member_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Multi", "Compara", "GeneMember");
my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Multi", "Compara", "Homology");
my $gene_member_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Multi", "Compara", "GeneMember");

print "Enter chromosome of variant: ";
my $chr_input = <STDIN>;
chomp($chr_input);

print "Enter bp location of variant in GR38: "; #based on version of API
my $bp_input = <STDIN>;
chomp($bp_input);

my $slice1 = $slice_adaptor->fetch_by_region("chromosome", $chr_input, $bp_input, $bp_input);
my $all_genes = $slice1->get_all_Genes;

foreach $gene ( @{$all_genes} ) {
    print "Variant location within " . $gene->external_name . " " . $gene->display_id . "\n";
    $base_Gid = $gene->display_id;
    $basegene = $gene->external_name;
}

my $gene_member = $gene_member_adaptor->fetch_by_stable_id($base_Gid);
my $homologies = $homology_adaptor->fetch_all_by_Member($gene_member);

foreach my $homology (@{$homologies}) {
	
	
	my @members = (@{$homology->get_all_Members});
	my $para_gene;	
	
	#Define hashes
	my %ENSPid = ();
	my %ENSTid = ();
	my %genename = ();
	my %geneobj = ();
	my %transslice = (); 
	my %strand = (); 
	my %fullseq = (); 
	my %trmapper = (); 
	my %peptide = ();
	my %peptide_start = ();

	
	#Are both members human?	
	my $hgmembercount = 0;
	foreach my $member (@members) {
		if ($member->taxon_id eq "9606"){
			$hgmembercount = $hgmembercount + 1;
		}
	}

	#If both members are human 		
	if ($hgmembercount == 2){
		
		my $simplealign = $homology->get_SimpleAlign(); 
		print $homology->description,": ";
		
		foreach my $example (@members) {

			my $ENSP = $example->stable_id;
			my $gene = $hg_adaptor->fetch_by_translation_stable_id($ENSP);
			if ($gene->external_name eq $basegene){
				$ENSPid{$basegene}=$example->stable_id; #confirm using longest trans for all
				$geneobj{$basegene}=$hg_adaptor->fetch_by_translation_stable_id($ENSP);
				$genename{$basegene}=$gene->external_name;
				$trans{$basegene}=$example->get_Transcript;
				$ENSTid{$basegene} = $trans{$basegene}->display_id;
			} else {
				$para_gene=$gene->external_name;
				$ENSPid{$para_gene}=$example->stable_id;
				$geneobj{$para_gene}=$hg_adaptor->fetch_by_translation_stable_id($ENSP);
				$genename{$para_gene}=$gene->external_name;
				$trans{$para_gene}=$example->get_Transcript;
				$ENSTid{$para_gene} = $trans{$para_gene}->display_id;			
			}
	
		}
		print "$genename{$para_gene} \n";

		$transslice{$basegene} = $slice_adaptor->fetch_by_transcript_stable_id($ENSTid{$basegene}); 
		$strand{$basegene} = $transslice{$basegene}->strand; 
		
		$trmapper{$basegene} = Bio::EnsEMBL::TranscriptMapper->new($trans{$basegene});   
		$trmapper{$para_gene} = Bio::EnsEMBL::TranscriptMapper->new($trans{$para_gene});  
		
		$fullseq{$basegene} = $simplealign->get_seq_by_id($ENSPid{$basegene});
		$fullseq{$para_gene} = $simplealign->get_seq_by_id($ENSPid{$para_gene});

		my %peptide_coord; 
		my @coordlist = $trmapper{$basegene}->genomic2pep($bp_input, $bp_input, $strand{$basegene}); 
		foreach my $var (@coordlist){
			$peptide_coord{$basegene} = $var; 
		}
	
		$peptide{$basegene} = $peptide_coord{$basegene}->start;
		my $col = $simplealign->column_from_residue_number($ENSPid{$basegene}, $peptide{$basegene});

		$peptide_coord{$para_gene} = $fullseq{$para_gene}->location_from_column($col);
		$peptide{$para_gene} = $peptide_coord{$para_gene}->start;

		if (defined $peptide{$para_gene}) { #how to deal with undefined?
			# print "$ENSPid{$basegene}-$peptide{$basegene} ~ $ENSPid{$para_gene}-$peptide{$para_gene}\n";
		} else {
			# print "$ENSPid{$basegene}-$peptide{$basegene} ~ Position - in $ENSPid{$para_gene}\n";
		}

		my %codoncoords;
		my @pos2bp_coords = $trmapper{$para_gene}->pep2genomic($peptide{$para_gene},$peptide{$para_gene});
		foreach my $var (@pos2bp_coords){
			$codoncoords{$basegene} = $var;
		}
		
		foreach my $seq ($simplealign->each_seq) {
			my $seqatlocation = $seq->subseq($col, $col);
			if ($seq->display_id eq $ENSPid{$basegene}) {
				print $genename{$basegene}, "\t", $ENSPid{$basegene}, "-", $peptide{$basegene}, ": ", $seqatlocation, "\n";
			} else {
				print $genename{$para_gene}, "\t", $ENSPid{$para_gene}, "-", $peptide{$para_gene}, ": ", $seqatlocation, "\n"
			}
		}
		
		print "\n";

		#***********************************# add filters below
		
		my $codon_start = $codoncoords{$basegene}->start;
		my $codon_end = $codoncoords{$basegene}->end;

		$transslice{$para_gene} = $slice_adaptor->fetch_by_transcript_stable_id($ENSTid{$para_gene});
		my $coord_sys2  = $transslice{$para_gene}->coord_system()->name();
		my $slice2_chr = $transslice{$para_gene}->seq_region_name();



		my $codon_slice2 = $slice_adaptor->fetch_by_region('chromosome', $slice2_chr, $codon_start, $codon_end);
		foreach my $vf ( @{ $vf_adaptor->fetch_all_by_Slice($codon_slice2) } ) {
			print $slice2_chr, ":", $vf->seq_region_start(), ' ', $vf->allele_string(), ' ', $vf->variation_name(), ' ', $vf->display_consequence, ' ';
			my $id = $vf->variation_name();
			my @csstates = @{$vf->get_all_clinical_significance_states};
			foreach my $sig ( @csstates ) {
				if ($sig =~ m/pathogenic/) {
					print "SIG: " , $sig, ' ';
				}
			}
		}
	
		print "\n\n************************\n"
	}
}	
		
		