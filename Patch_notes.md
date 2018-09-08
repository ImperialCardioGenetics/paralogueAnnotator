Patch notes:

7/8/18 
Added: 
```
Bio::EnsEMBL::Registry->set_reconnect_when_lost(1);
```
to ensure connection to Compara is not lost due to inactivity.

6/9/18
Added:
```
next if (!defined $col);

next if (!defined $fullseq{$para_gene});

next if (!defined $peptide_coord{$para_gene});
```
to error catch "Can't call start/location_from_column methods" error.

7/9/18
Added:
```
$genemember_adaptor->dbc->disconnect_if_idle;

$hg_adaptor->dbc->disconnect_if_idle;

$slice_adaptor->dbc->disconnect_if_idle;

$variationfeature_adaptor->dbc->disconnect_if_idle;

$transcriptvariation_adaptor->dbc->disconnect_if_idle;

$homology_adaptor->dbc->disconnect_if_idle;
```
to prevent "MySQL server has gone away" error