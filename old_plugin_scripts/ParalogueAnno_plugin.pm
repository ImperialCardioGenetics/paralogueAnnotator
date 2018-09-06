package ParalogueAnno_plugin;
#APIs installed this way will automatically be on the release/78 branch.
#how to look at only missense?
#test that all variables are actually used

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::LocatableSeq;
use Bio::EnsEMBL::TranscriptMapper;	


use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);


sub feature_types {
        return ['Transcript'];  #yes, use transcript since we only want missense (aka trans altering) http://www.ensembl.org/info/docs/Doxygen/variation-api/classBio_1_1EnsEMBL_1_1Variation_1_1Utils_1_1BaseVepPlugin.html#a1c1872ebd38c916ac0850c870a48bdbe
    }

sub variant_feature_types {
    return ['VariationFeature'];
}

sub get_header_info {
	return {
		Paralogue_Vars => "List of pathogenic variants in paralogue genes"
	};
}

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);

    my $config = $self->{config};
    my $reg = 'Bio::EnsEMBL::Registry';

    # reconnect to DB without species param
    if($config->{host}) {
        $reg->load_registry_from_db(
            -host       => $config->{host},
            -user       => $config->{user},
            -pass       => $config->{password},
            -port       => $config->{port},
            -db_version => $config->{db_version},
            -no_cache   => $config->{no_slice_cache},
        );
    }
	#$db = $self->{config}->{db_version};
	
    $self->{homology_adap} = $config->{homo} ||= $reg->get_adaptor('Multi', 'compara', 'Homology')
        or die "Failed to fetch conservation adaptor\n";

    return $self;
}


sub run {

	my ($self, $tva) = @_;
	my $btv = $tva->base_transcript_variation();
	my $vf = $btv->variation_feature;
	
	my $lineid = ${$btv->transcript}{"stable_id"}; #### if this field (trans id) is equiv to homology trans, the return result
	my $result = "";
	
	if ($vf->display_consequence() eq 'missense_variant' && $lineid) {
	
		my $genome_db_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'GenomeDB');
		my $hg_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human","Core","Gene"); 
		my $familyA = Bio::EnsEMBL::Registry->get_adaptor("Multi", "Compara", "Family");
		my $seqmemberA = Bio::EnsEMBL::Registry->get_adaptor("Multi", "Compara", "SeqMember");
		my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human", "Core", "Slice");
		my $transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human", "Core", "Transcript");
		my $vf_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human", "Variation", "Variationfeature");
		my $pa = Bio::EnsEMBL::Registry->get_adaptor("Human","Variation","Phenotype");
		my $gene_member_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Multi", "compara", "GeneMember");
		my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');

		my @variants = ();
		my %variation_gene = ();
		my %variation_chr = ();
		my %variation_start = ();
		my %variation_allele = ();
		my %variation_name = ();
		my %variation_cons = ();
		my $base_Gid = ();
		my $basegene = ();
		my $para_gene = (); 	
		
		my $col = ();
	
	
		my $chr_input = $vf->seq_region_name();
		my $bp_input = $vf->start();

		my $slice1 = $slice_adaptor->fetch_by_region("chromosome", $chr_input, $bp_input, $bp_input);
		my $all_genes = $slice1->get_all_Genes;

		foreach my $example ( @{$all_genes} ) { #find a different way to do this
			$base_Gid = $example->display_id;
			$basegene = $example->external_name;
#			$listofgenes .= $base_Gid;
		}
		
 		my $gene_member = $gene_member_adaptor->fetch_by_stable_id($base_Gid);
 		my $homologies = $self->{homology_adap}->fetch_all_by_Member($gene_member);

		my @homologyArray = (@{$homologies});		
#	  		foreach my $homology (@homologyArray) {
		foreach my $homology (@{$homologies}) {
	
			my @members = (@{$homology->get_all_Members});

				#Define hashes
			my %ENSPid = ();
			my %ENSTid = ();
			my %genename = ();
			my %geneobj = ();
			my %transslice = (); 
			my %strand = (); 
			my %fullseq = (); 
			my %trmapper = (); 
			my %peptide_coord = (); 
			my %peptide = ();
			my %peptide_start = ();
			my %resatlocation = ();
			my %trans = ();
	
		
				#Are both members human?	
			my $hgmembercount = 0;
			foreach my $member (@members) {							
				if ($member->taxon_id eq "9606"){
					$hgmembercount = $hgmembercount + 1;
				}
			}
		
			if ($hgmembercount == 2){ 	
				foreach my $member (@members) {			
					my $ENSP = $member->stable_id; #confirm using longest trans for all
					my $gene = $hg_adaptor->fetch_by_translation_stable_id($ENSP);					
					if ($gene->external_name eq $basegene){
						$ENSPid{$basegene}=$member->stable_id; #confirm using longest trans for all
						$geneobj{$basegene}=$hg_adaptor->fetch_by_translation_stable_id($ENSP);
						$genename{$basegene}=$gene->external_name;
						$trans{$basegene}=$member->get_Transcript;
						$strand{$basegene}=$trans{$basegene}->strand;
						$ENSTid{$basegene} = $trans{$basegene}->display_id;
						$trmapper{$basegene} = Bio::EnsEMBL::TranscriptMapper->new($trans{$basegene});  
					} else {
						$para_gene=$gene->external_name;
						$ENSPid{$para_gene}=$member->stable_id;
						$geneobj{$para_gene}=$hg_adaptor->fetch_by_translation_stable_id($ENSP);
						$genename{$para_gene}=$gene->external_name;
						$trans{$para_gene}=$member->get_Transcript;
						$ENSTid{$para_gene} = $trans{$para_gene}->display_id;	
						$trmapper{$para_gene} = Bio::EnsEMBL::TranscriptMapper->new($trans{$para_gene});
					}								
				}
				if ($lineid eq $ENSTid{$basegene}) { 
					my $simplealign = $homology->get_SimpleAlign();
					$fullseq{$basegene} = $simplealign->get_seq_by_id($ENSPid{$basegene});
					$fullseq{$para_gene} = $simplealign->get_seq_by_id($ENSPid{$para_gene});	
					my @coordlist = $trmapper{$basegene}->genomic2pep($bp_input, $bp_input, $strand{$basegene}); #when list has one element how to extract?
 					foreach my $var (@coordlist){ 
 						$peptide_coord{$basegene} = $var; 
 					}
 					$peptide{$basegene} = $peptide_coord{$basegene}->start;
 					$col = $simplealign->column_from_residue_number($ENSPid{$basegene}, $peptide{$basegene});
 					$peptide_coord{$para_gene} = $fullseq{$para_gene}->location_from_column($col);
 					$peptide{$para_gene} = $peptide_coord{$para_gene}->start; 
  					my %codoncoords;
 					my @pos2bp_coords = $trmapper{$para_gene}->pep2genomic($peptide{$para_gene},$peptide{$para_gene});
 					foreach my $var (@pos2bp_coords){
 						$codoncoords{$para_gene} = $var;
 					}
 					my $codon_start = $codoncoords{$para_gene}->start;
 					my $codon_end = $codoncoords{$para_gene}->end; 		
 					
 					$transslice{$para_gene} = $slice_adaptor->fetch_by_transcript_stable_id($ENSTid{$para_gene});
 					my $coord_sys2  = $transslice{$para_gene}->coord_system()->name();
 					my $slice2_chr = $transslice{$para_gene}->seq_region_name(); 					
 					my $codon_slice2 = $slice_adaptor->fetch_by_region('chromosome', $slice2_chr, $codon_start, $codon_end);
 					
 					foreach my $vf ( @{ $vf_adaptor->fetch_all_by_Slice($codon_slice2) } ) {
 						my @csstates = @{$vf->get_all_clinical_significance_states};
						my $path = 0;
						foreach my $sig ( @csstates ) {
							if ($sig =~ m/pathogenic/) {
									$path++;
							}
						}
						if ($path > 0 || $vf->allele_string() eq "HGMD_MUTATION"){ 					
							my $id = $vf->variation_name();
							push @variants, $id;
							$variation_gene{$id} = $genename{$para_gene};
							$variation_name{$id} = $vf->variation_name();
						}							 						
 					}
																																					
				}		
			}
 		}
		foreach my $ids (@variants) {
			$result .= $variation_gene{$ids} . " " . $variation_name{$ids} . "; "; 
		} 	
 	}			
		
	return {
	#	Paralogue_Vars => $results
		Paralogue_Vars => $result
	};				
			
} 

1;
