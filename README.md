
# Pathway Tool -- Tool for Pathway enrichment analysis

## Donwload and Installation

Download pathway tool module from [github](https://github.com/tripal/tripal_pathway)
```
cd /var/www/html/youSiteFolder/sites/all/modules/
git clone https://github.com/tripal/tripal_pathway
```

Install pathway tool through "Administration of Modules" of Drupal.

Or install pathway tool using command:
```
drush pm-enable tripal_pathway
```

## Add Pathway File

The pathway file is produced by processing the output of [PathwayTools](http://brg.ai.sri.com/ptools/). 
First, we use PathwayTools to predict the pathways for the annotated genes. The will generate a file
named pathways.col includes predict pathways for genes. Next, we retrieve information from pathways.col
to generate the pathway file. The pathway file is tab-delimit file which include 4 columns. 
The format of pathway file:
- Gene ID
- Functional Description
- Pathway ID
- Enzyme Name

Please check the example of pathway file under: **example/wm_97103_kept.pwy**

Then load the pathway file to database.
> Home -> Add Content -> Pathway File  


