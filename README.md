# paralogueAnnotator

The repository is the home of an ensembl VEP plug in for genomic DNA variation.  Given details of a variant, the plug-in determines whether variation has been previously observed at an equivalent position in other members of the protein family.  If so, details of the variants in paralogues are returned, including phenotypic annotations.  This approach has been shown to identify human disease-associated variants across protein families.

See http://cardiodb.org/Paralogue_Annotation/ for an overview of the approach, with specific applications to inherited cardiac conditions.

Previous publications:  

- [Paralogous annotation of disease-causing variants in Long QT syndrome genes](http://onlinelibrary.wiley.com/doi/10.1002/humu.22114/abstract)  
- [Paralogue annotation identifies novel pathogenic variants in patients with Brugada syndrome and catecholaminergic polymorphic ventricular tachycardia](http://jmg.bmj.com/content/early/2013/10/17/jmedgenet-2013-101917.full)  

Example and recommended usage:

General command line usage are as follows. For optimal performance, we suggest running VEP offline with a local Ensembl cache. 
For `variants` mode:
```
perl -I [directory of plugin] [directory of installed VEP] --force_overwrite --vcf --offline --cache --dir_cache [directory of cache] -i [input file path] -o [output file path] --plugin ParalogueAnno_plugin_cleanup
```
For `paraloc` mode:
```
perl -I [directory of plugin] [directory of installed VEP] --force_overwrite --vcf --offline --cache --dir_cache [directory of cache] -i [input file path] -o [output file path] --plugin ParalogueAnno_plugin_cleanup,paraloc
```

_VEP version 90 used in example below_


```
perl -I /home/user /data/Install/ensembl-vep/vep --force_overwrite --vcf --offline --cache --dir_cache /home/user -i /home/user/input_file.vcf -o /home/user/output_file.out --plugin ParalogueAnno_plugin_cleanup
```

For documentation on getting VEP to run faster please visit http://www.ensembl.org/info/docs/tools/vep/script/vep_other.html