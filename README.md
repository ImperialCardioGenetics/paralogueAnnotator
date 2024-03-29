# paralogueAnnotator

The repository is the home of an ensembl VEP plug in for genomic DNA variation.  Given details of a variant, the plug-in determines whether variation has been previously observed at an equivalent position in other members of the protein family.  If so, details of the variants in paralogues are returned, including phenotypic annotations.  This approach has been shown to identify human disease-associated variants across protein families.

See http://cardiodb.org/Paralogue_Annotation/ for an overview of the approach, with specific applications to inherited cardiac conditions.

Previous publications:  

- [Paralogous annotation of disease-causing variants in Long QT syndrome genes](http://onlinelibrary.wiley.com/doi/10.1002/humu.22114/abstract)  
- [Paralogue annotation identifies novel pathogenic variants in patients with Brugada syndrome and catecholaminergic polymorphic ventricular tachycardia](http://jmg.bmj.com/content/early/2013/10/17/jmedgenet-2013-101917.full)  

### Installation and setup:

Install Ensembl's [VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) accordingly.

**Trouble Shooting notes:**
Perl modules `DBI` and `DBD::mysql` may not be installed in recognised directory paths when using `cpanm` if there are multiple perl versions on the system you're using. We recommend using [Perlbrew](https://perlbrew.pl/) to manage and switch between different versions of perl. Then try installing modules directly via:
```
perl -MCPAN -e 'install DBD::mysql'
perl -MCPAN -e 'install DBI'
```

Make sure the perl and perl API environments are setup correctly before running, especially on systems that do not permanently save user profiles. For example, this can be done by either typing into bash command line the following:
```
PERL5LIB=${PERL5LIB}:/home/user/BioPerl-1.6.924
PERL5LIB=${PERL5LIB}:/home/user/ensembl/modules
PERL5LIB=${PERL5LIB}:/home/user/ensembl-compara/modules
PERL5LIB=${PERL5LIB}:/home/user/ensembl-variation/modules
PERL5LIB=${PERL5LIB}:/home/user/ensembl-funcgen/modules
PERL5LIB=${PERL5LIB}:/home/user/.cpanm
PERL5LIB=${PERL5LIB}:/home/user/.cpan/build
export PERL5LIB
```
Or by creating a `.bashrc` file in your home directory with the above and running `source .bashrc`. Change `/home/user/` to any appropriate installed directories if required.


### How to run VEP with perl plugin:

General command line usage are as follows. 
```
perl -I [directory of plugin] [directory of installed VEP] \
    --force_overwrite --vcf --offline --cache --dir_cache [directory of cache] \
    -i [input file path] -o [output file path] \
    --plugin ParalogueAnnotation --assembly [assebmly build] --port [port number]
```

### What the arguments indicate:

`-I`: needed to point to location of the plugin (adding to @INC)

`--force_overwrite`: allows you to overwrite the output file

`--vcf`: the format of the input file

`--offline` & `--cache`: both attempt to use the cache, although the offline mode will override when it can't access Compara data

`--dir_cache`: location of our data cache

`-i`: input file

`-o`: output file

`--plugin`: use name of plugin, enter plugin inputs separated by commas (see below for possible options)

`--assembly`: indicates the assembly build of the genomic coordinates of input vcf. Use either `GRCh38` or `GRCh37`. If not specified then `GRCh38` is used. If using `GRCh37` then must be used with `--port 3337` (see below).

`--port`: use `--port 3337` if using GRCh37 genome coordinates in order for compara to also use GRCh37 coordinates. If no port is specified then default GRCh38 is used.

For optimal performance, we suggest running VEP offline with a local Ensembl cache. For documentation on getting VEP to run faster please visit http://www.ensembl.org/info/docs/tools/vep/script/vep_other.html

### Additional `--plugin` options:

`[--plugin ParalogueAnnotation]`(_default_) or `[--plugin ParalogueAnnotation,variants,damaging]`: VEP plugin will run in __variant__ mode, outputting damaging variants that appear in every paralogue's equivalent location

`[--plugin ParalogueAnnotation,variants,all]`: VEP plugin will run in __variant__ mode, outputting ALL variants that appear in every paralogue's equivalent location

`[--plugin ParalogueAnnotation,paraloc]`: VEP plugin will run in __location__ mode, outputting the variant's equivalent location (3bp represents AA codon) in each paralogue

For example:
```
perl -I /home/user /data/Install/ensembl-vep/vep --force_overwrite --vcf --offline --cache --dir_cache /home/user -i /home/user/input_file.vcf -o /home/user/output_file.out --plugin ParalogueAnnotation --assembly GRCh38
```
(_VEP version 90 used in example above_)


In __location__ mode (option=`paraloc`) the Paralogue_Vars info field gives the following information:

    ParaSym: Symbol of Paralogue Gene
    ParaChr: Chromosome of Paralogue Gene
    ParaStart: Start location of equivalent AA codon in paralogue gene
    ParaEnd: End location of equivalent AA codon in paralogue gene
    AAorig: Reference amino acid in original gene
    AApara: Reference amino acid in paralogue gene
    REFID: Binary indicating whether both reference amino acids are the same

    in the format:
    |ParaSym:ParaChr_ParaStart-ParaEnd:AAorig:AApara:REFID|

    example of "paraloc" output for variant chr3:38630420 (in SCN5A)
    |SCN10A:chr3_38792157-38792159:V:V:REFID=1|

In __variant__ mode (option=`variants`) the Paralogue_Vars info field gives the following information:

    ParaSym: Symbol of Paralogue Gene
    ParaChr: Chromosome of Paralogue Gene
    ParaStart: Start location of equivalent AA codon in paralogue gene
    ParaEnd: End location of equivalent AA codon in paralogue gene
    VarID: ID of the variant in paralogue
    REFID: Binary indicating whether both reference amino acids are the same
    origAAs: Amino acid string for original variant as "REF/ALT"
    paraAAs: Amino acid string for paralogue variant as "REF/ALT"

    in the format:
    |ParaSym:ParaChr_ParaStart-ParaEnd:VarID:REFID:origAAs:paraAAs|


    example of "variants" output for variant chr3:38630420 (in SCN5A)
    |SCN10A:chr3_38792158:rs202143516:REFID=1:V/I:V/G|

In __location__ mode there will be only one output per paralogue gene; in __variant__ mode there may by multiple variants (and therefore outputs) per gene.

Paralogue variants that are output from this plugin come from dbSNP (which includes ClinVar variants) and HGMD: variants beginning with "rs" come from dbSNP; variants beginning with "C" come from HGMD.

The `damaging` option shows variants that have *at least* one annotation as "damaging" in the ensembl database.
