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
my $hg_adaptor = $reg->get_adaptor("human","core","Gene"); 
my $familyA = $reg->get_adaptor("Multi", "compara", "Family");
my $seqmemberA = $reg->get_adaptor("Multi", "compara", "SeqMember");
my $slice_adaptor = $reg->get_adaptor("Human", "core", "Slice");
my $transcript_adaptor = $reg->get_adaptor("Human", "core", "Transcript");
my $vf_adaptor = $reg->get_adaptor('human', 'variation', 'variationfeature');
my $pa = $reg->get_adaptor("human","variation","phenotype");


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
    $seq1Pid = $gene->display_id;
}

my $gene_member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'GeneMember');
my $gene_member = $gene_member_adaptor->fetch_by_stable_id($seq1Pid);


my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');
my $homologies = $homology_adaptor->fetch_all_by_Member($gene_member);

my $n = 0;

foreach my $homology (@{$homologies}) {
		
		my $count = 0;
		my $counthuman = 0;
		my @members = (@{$homology->get_all_Members});
		
		foreach my $example (@members) {
			$count = $count+1;
			my $taxon = $example->taxon_id;
				if ($taxon eq "9606") {
					$counthuman = $counthuman+1;
				}
		}

		if ($count==$counthuman){
			print $n," ", $homology->description,"\n";
			foreach my $example (@members) {
				my $ENSP = $example->stable_id; #confirm using longest trans for all
				my $gene = $hg_adaptor->fetch_by_translation_stable_id($ENSP);
				my $trans = $example->get_Transcript;
				print "\t" . $gene->external_name . " " . $example->stable_id . "\n";
			}
		}
			
		$count = 0;
		$counthuman = 0;

$n = $n+1;
}

print "Which homology (enter number) \n"; 
my $homology_input = <STDIN>; 
chomp($homology_input); 

my $selectedhomology = $homologies->[$homology_input]; 
my $simplealign = $selectedhomology->get_SimpleAlign();
my $simplealignCDS = $selectedhomology->get_SimpleAlign(-seq_type => 'cds');

my @selmembers = (@{$selectedhomology->get_all_Members});

foreach my $example (@selmembers) {
				my $ENSP = $example->stable_id;
				my $gene = $hg_adaptor->fetch_by_translation_stable_id($ENSP);
				my $trans = $example->get_Transcript;
				print "\t" . $gene->external_name . " " . $example->stable_id . "\n";
}

my $alncount = 0;
my $seq1Pid;
my $seq2Pid;
my $gene1;
my $gene2;
my $trans1;
my $trans2;
my $seq1Tid;
my $member1Pid;
my $member2Pid;
my $member1Tid;
my $member2Tid;

my $membercount = 0;

my @selmembers = (@{$selectedhomology->get_all_Members});
foreach my $selmember (@selmembers) {
	$membercount=$membercount+1;
	if ($membercount==1){
		$member1Pid = $selmember->stable_id;
		$gene1 = $hg_adaptor->fetch_by_translation_stable_id($member1Pid);
		$trans1 = $selmember->get_Transcript;
		$member1Tid = $trans1->display_id;
		$trans1b = $transcript_adaptor->fetch_by_stable_id($member1Tid);
	}
	if ($membercount==2){
		$trans2 = $selmember->get_Transcript;
		$member2Tid = $trans2->display_id;
		$member2Pid = $selmember->stable_id;
		$gene2 = $hg_adaptor->fetch_by_translation_stable_id($member2Pid);
	}
}

my $slice1 = $slice_adaptor->fetch_by_transcript_stable_id($member1Tid);
my $strand1 = $slice1->strand;

print $gene1->external_name . " transcript id: $member1Tid, peptide id: $member1Pid, strand: $strand1\n";
print $gene2->external_name . " transcript id: $member2Tid, peptide id: $member2Pid \n";

foreach $seq ( $simplealignCDS->each_seq() ) {
	$alncount=$alncount+1;
	if ($alncount==1){
		$seq1Pid=$seq->id;
	}
	if ($alncount==2){
		$seq2Pid=$seq->id;
	}
}

if ($seq1Pid eq $member1Pid && $seq2Pid eq $member2Pid){
	print "Alignment order matches member order.\n\n";
} 


my $fullseq1 = $simplealign->get_seq_by_id($seq1Pid);
my $fullseq2 = $simplealign->get_seq_by_id($seq2Pid);

my $trmapper1 = Bio::EnsEMBL::TranscriptMapper->new($trans1);
my $trmapper2 = Bio::EnsEMBL::TranscriptMapper->new($trans2);

my $coordhash1;

my @coordlist1 = $trmapper1->genomic2pep($bp_input, $bp_input, $strand1); 
	foreach my $var1 (@coordlist1){
		$coordhash1 = $var1;
	}

my $pos1pep = $coordhash1->start;

my $col = $simplealign->column_from_residue_number($seq1Pid, $pos1pep);

$pos2pep = $fullseq2->location_from_column($col);

if (defined $pos2pep){
	print "Peptide position $pos1pep in $seq1Pid : Column $col of alignment : Peptide position " . $pos2pep->start . " in $seq2Pid\n";
} else {
	print "Position $pos1pep in $seq1Pid : Column $col of alignment : Position - in $seq2Pid\n";
}

my $pos2pep_start = $pos2pep->start;


my $coordhash;
my @pos2bp_coords = $trmapper2->pep2genomic($pos2pep_start,$pos2pep_start);
foreach my $var (@pos2bp_coords){
	$coordhash = $var;
	#print $coordhash;



}

my $slice2 = $slice_adaptor->fetch_by_transcript_stable_id($member2Tid);
my $coord_sys2  = $slice2->coord_system()->name();
my $slice2_chr = $slice2->seq_region_name();

my $codon_start = $coordhash->start;
my $codon_end = $coordhash->end;
	
print "Genomic coords of homologous codon in $seq2Pid: " . $coord_sys2 . " " . $slice2_chr . ":" . $coordhash->start . "-" . $coordhash->end . "\n";	
	
my %count;
foreach my $seq ($simplealign->each_seq) {
   	my $seqatlocation = $seq->subseq($col, $col);
   	$count{$seqatlocation}++;
   	my $name = $seq->display_id;
   	print $name,": ", $seqatlocation, "\n";
}

#foreach $res (keys %count) {
#	printf "$res: $count{$res}\n";
#}	

print "\n";

my $codon_slice2 = $slice_adaptor->fetch_by_region('chromosome', $slice2_chr, $codon_start, $codon_end);
foreach my $vf ( @{ $vf_adaptor->fetch_all_by_Slice($codon_slice2) } ) {
    print $slice2_chr, ":", $vf->seq_region_start(), ' ', $vf->allele_string(), ' ', $vf->variation_name(), ' ', $vf->display_consequence, ' ';
	my $id = $vf->variation_name();
	my @csstates = @{$vf->get_all_clinical_significance_states};
	foreach my $sig ( @csstates ) {
		print "SIG: " , $sig, ' ';
	}
	print "\n";
#		my @csstates = @{$vf->get_all_clinical_significance_states};
#		foreach my $sig ( @csstates ) {
#			print $sig, ' ';
#		}	
	
 }
 print "\n";
 
 
#print $codon_slice2->clinical_significance;