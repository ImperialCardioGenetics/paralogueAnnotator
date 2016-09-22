#!/usr/bin/perl

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::AlignIO;
use Bio::LocatableSeq;

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


print "Enter Gene ID \n";
my $gene_input = <STDIN>;
chomp($gene_input);


my $gene_member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'GeneMember');
my $gene_member = $gene_member_adaptor->fetch_by_stable_id($gene_input);

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
				my $ENSP = $example->stable_id;
				my $gene = $hg_adaptor->fetch_by_translation_stable_id($ENSP);
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


my $alignIO = Bio::AlignIO->newFh(
        -interleaved => 0,
        -fh          => \*STDOUT,
        -format      => "clustalw",
        -idlength    => 20);

my $alncount = 0;
my $seq1;
my $seq2;


foreach $seq ( $simplealignCDS->each_seq() ) {
	$alncount=$alncount+1;
	if ($alncount==1){
		$seq1=$seq->id;
	}
	if ($alncount==2){
		$seq2=$seq->id;
	}
}

$fullseq1 = $simplealignCDS->get_seq_by_id($seq1);
$fullseq2 = $simplealignCDS->get_seq_by_id($seq2);

print "Which bp in $seq1 cds? (enter number): ";

while(my $pos1 = <>){
	chomp($pos1);
	$col = $simplealignCDS->column_from_residue_number($seq1, $pos1);

	$pos2 = $fullseq2->location_from_column($col);
	if (defined $pos2){
	print "Position $pos1 in $seq1 : Column $col of alignment : Position " . $pos2->start . " in $seq2\n";
	} else {
	print "Position $pos1 in $seq1 : Column $col of alignment : Position - in $seq2\n";
	}
	
	my %count;
	foreach my $seq ($simplealignCDS->each_seq) {
    	my $seqatlocation = $seq->subseq($col, $col);
    	$count{$seqatlocation}++;
    }
    foreach $res (keys %count) {
		printf "$res: $count{$res}\n";
	}
	
	print "\nWhich bp in $seq1 cds? (enter number): ";
}