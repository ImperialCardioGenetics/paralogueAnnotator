package ParalogueAnno_plugin_cleanup;
#APIs installed this way will automatically be on the release/78 branch.
#how to look at only missense?
#test that all variables are actually used
#	first parameter, "paraloc" for paralog variant location, or "variants"; variants is default
#   second parameter, "all" for all variants or "damaging" for only damaging variant, "vcf" to search file; damaging is default
#   third input vcf file to search


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
		Paralogue_Vars => "Equivalant variants and locations in paralogous genes"
	};
}

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);

    my $params = $self->params;
   
    shift @$params if $params->[0] && $params->[0] eq '1';  # REST API passes 1 as first param

	$self->{run} = $params->[0] || 'variants';
    $self->{output} = $params->[1] || 'damaging';
    $self->{file} = $params->[2];

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
	
    $self->{homology_adap} = $config->{homo} ||= $reg->get_adaptor('Multi', 'compara', 'Homology')
        or die "Failed to fetch conservation adaptor\n";

    return $self;
}


sub run {

	my ($self, $tva) = @_;
	my $result = "";
	
	my $btv = $tva->base_transcript_variation();
	my $tv = $tva->transcript_variation;
	my $vf = $btv->variation_feature;
	
	my $basepepchange = $tv->pep_allele_string;
	my $lineid = ${$btv->transcript}{"stable_id"}; 
	
	
	if ($vf->display_consequence() eq 'missense_variant'){# && $lineid) {


		#Define adaptors
	
		my $genome_db_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'GenomeDB');
		my $hg_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human","Core","Gene"); 
		my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human", "Core", "Slice");
		my $transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human", "Core", "Transcript");
		my $variationfeature_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human", "Variation", "Variationfeature");
		my $transcriptvariation_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Human", "Variation", "TranscriptVariation");
		my $genemember_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Multi", "compara", "GeneMember");
		my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');

		#Define arrays and hashes

		my @variants = ();
		my %variation_gene = ();
		my %variation_chr = ();
		my %variation_start = ();
		my %variation_allele = ();
		my %variation_name = ();
		my %variation_cons = ();
		my %variation_pos = ();
		my %variation_pepchange = ();
		my %variation_pep = ();
		my %tvariation = ();
		my %variation_REFid = ();
		my %varREFbase = ();
		my %varREFpara = ();
		my %REFresatlocation = ();
		my %ALTresatlocation = ();
		
		#Define variables
		
		my $base_Gid = ();		
		my $basegene = ();
		my $para_gene = ();
		my $REFid = (); 
		my $col = ();
	
	
		#Define variant input
	
		my $chr_input = $vf->seq_region_name();
		my $bp_input = $vf->start();
		my $slice_input = $slice_adaptor->fetch_by_region("chromosome", $chr_input, $bp_input, $bp_input);
		my $all_genes = $slice_input->get_all_Genes;

		foreach my $example ( @{$all_genes} ) { #find a different way to do this
			$base_Gid = $example->display_id;
			$basegene = $example->external_name;
		}
		
		#Fetch all homologies
		
 		my $gene_member = $genemember_adaptor->fetch_by_stable_id($base_Gid);
 		my $homologies = $self->{homology_adap}->fetch_all_by_Member($gene_member);

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
			my %trans = ();
	
			#Are both genes human?	
			
			my $hgmembercount = 0;
			foreach my $member (@members) {							
				if ($member->taxon_id eq "9606"){
					$hgmembercount = $hgmembercount + 1;
				}
			}
			
			#For homologies with two human genes, set basegene and para_gene
			
			if ($hgmembercount == 2){ 	
				foreach my $member (@members) {			
					my $ENSP = $member->stable_id; #confirm using longest trans for all
					my $gene = $hg_adaptor->fetch_by_translation_stable_id($ENSP);					
					if ($gene->external_name eq $basegene){
						$ENSPid{$basegene}=$member->stable_id;
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
				
				#Find paralogous location using transcripts aligned in stored EnsEMBL homology
				
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
 					foreach my $var (@pos2bp_coords){ #extract one element
 						$codoncoords{$para_gene} = $var;
 					}
 					my $codon_start = $codoncoords{$para_gene}->start;
 					my $codon_end = $codoncoords{$para_gene}->end; 		
 					
 					$transslice{$para_gene} = $slice_adaptor->fetch_by_transcript_stable_id($ENSTid{$para_gene});
 					my $coord_sys2  = $transslice{$para_gene}->coord_system()->name();
 					my $slice2_chr = $transslice{$para_gene}->seq_region_name(); 					
 					my $codon_slice2 = $slice_adaptor->fetch_by_region('chromosome', $slice2_chr, $codon_start, $codon_end);

 					$REFresatlocation{$basegene} = $fullseq{$basegene}->subseq($col, $col);
 					$REFresatlocation{$para_gene} = $fullseq{$para_gene}->subseq($col, $col); 
 					my %REFid = ();
 					
 					if ($REFresatlocation{$basegene} eq $REFresatlocation{$para_gene}) {
						$REFid{$para_gene} = 1;
					} else {
						$REFid{$para_gene} = 0;
					}
 					
 					if ($self->{run} eq "paraloc") { 
 						$result .= "|$para_gene:chr$slice2_chr" . "_$codon_start-$codon_end:$REFresatlocation{$basegene}:$REFresatlocation{$para_gene}:REFID=$REFid{$para_gene}";					
 						next;	
 					}
 					
# 					if $self->{run} eq "vcf" {
# 						use Tabix;
# 						my $tabix = Tabix->new('-data' => '/data/Mirror/gnomAD/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.14.vcf.gz'); 
# 						my $var = $tabix->read(
#   						$tabix->query( 14, 23415078, 23415078) 
# 						); 				
# 						$result = $var;
# 						next;
# 					}		
					
#					 if $self->{run} eq "vcf" {
# 						use Tabix;
# 						my $tabix = Tabix->new('-data' => $file); 
# 						$var = $tabix->read(
#   						$tabix->query( $slice2_chr, $codon_start, $codon_end) 
# 						); 				
# 						$result = $var;
# 						next;
# 					}	
								
 					foreach my $vf ( @{ $variationfeature_adaptor->fetch_all_by_Slice($codon_slice2) } ) {
 						my @csstates = @{$vf->get_all_clinical_significance_states};
						my $path = 0;
						foreach my $sig ( @csstates ) {
							if ($sig =~ m/pathogenic/) {
									$path++;
							}
						}	
						
						my $id = ();
					
						my $passvars = sub {
							my @array = ();
							$array[0] = $_[0];
							$id = $_[0]->variation_name();
							push @variants, $id;
							$tvariation{$id} = $transcriptvariation_adaptor->fetch_all_by_VariationFeatures(\@array)->[0];
							$variation_gene{$id} = $genename{$para_gene};
							$variation_name{$id} = $_[0]->variation_name();
							$variation_cons{$id} = $_[0]->display_consequence();
							$variation_chr{$id} = $slice2_chr;
							$variation_pos{$id} = $_[0]->seq_region_start();
							$variation_pepchange{$id} = $tvariation{$id}->pep_allele_string();
							$variation_pep{$id} = (split $variation_pepchange{$id},"/")[1];
							$varREFbase{$id} = $REFresatlocation{$basegene};
							$varREFpara{$id} = $REFresatlocation{$para_gene};
							if ($REFresatlocation{$basegene} eq $REFresatlocation{$para_gene}) {
								$variation_REFid{$id} = 1;
							} else {
								$variation_REFid{$id} = 0;
							}
						}; 					
						
						if ($vf->var_class eq "SNP") {
						
							if ($self->{output} eq "all") {
								&$passvars($vf);
							} elsif ($self->{output} eq "damaging" && ($path > 0 || $vf->allele_string() eq "HGMD_MUTATION")){
								&$passvars($vf);
							}											 						
 						}  					
 					}	 																																									
				}
			}		
 		}
		
		foreach my $ids (@variants) {
			$result .= "|$variation_gene{$ids}:chr$variation_chr{$ids}_$variation_pos{$ids}:$variation_name{$ids}:REFID=$variation_REFid{$ids}:$basepepchange:$variation_pepchange{$ids}"; 
		}
		
		if ($result ne "") {		 			
			return {
				Paralogue_Vars => $result . "|"
			};				
		} else {		#so no errors are thrown at runtime for empty results hash
			return {
			};
		}				
 	} else {		#so no errors are thrown at runtime for empty results hash
		return {
		};
 	}	
 	
 			
} 

1;
