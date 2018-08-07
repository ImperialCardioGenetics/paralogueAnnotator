Patch notes:

7/8/18 - Added line 14: 
```
Bio::EnsEMBL::Registry->set_reconnect_when_lost();
```
to ensure connection to Compara is not lost due to inactivity.