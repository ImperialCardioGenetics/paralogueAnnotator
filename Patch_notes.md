Patch notes:

7/8/18 
Added line 14: 
```
Bio::EnsEMBL::Registry->set_reconnect_when_lost(1);
```
to ensure connection to Compara is not lost due to inactivity.

6/9/18
Added line 227:
```
next if (!defined $col);
```
Added line 228:
```
next if (!defined $fullseq{$para_gene});
```
Added line 230:
```
next if (!defined $peptide_coord{$para_gene});
```

7/9/18
Added line 98:
```
$genemember_adaptor->dbc->disconnect_if_idle;
```