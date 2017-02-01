<?php
/**
 * This template displays the help page for the go tool
 */
?>

<h3>Module Description</h3>
<p>This module provides a basic interface to allow your users to perform pathway enrichment analysis.</p>

<h3>Setup Instructions</h3>
<ol>
<li>Prepare pathway annotation file </li>
The pathway annotation file is a tab delimit file with 4 columns: <br>
<code>gene_id, gene_description(option), pathway_id, pathway_name </code>

<li>Install perl module Statistics::Multtest</li>
<code>sudo perl -MCPAN -e shell <br>
  install Statistics::Multtest
</code>

<li>create soft link for pwy_tool.pl</li>
Example: 
<code>ln -s /var/www/html/sites/all/modules/tripal_pathway/pwy_tool.pl /usr/local/bin</code>

<li><a href="<?php print url('node/add/pwy');?>">Create "Pathway File" nodes</a> for each dataset</li>
	<ul>
	<li>Human-readable name for pathway annotation file. For example, watermelon genome (97103) means the pathway annotation file is predicted using all genes of watermelon genome (97103).</li>
	<li>Locateion of pathway annotation file, please use full path here. Example: /var/www/html/sites/all/modules/tripal_pathway/example/wm_97103_kept.pwy</li>
    <li>The url link of gene ID. This link should not includes "http://" on the left, and must include back slash "/" if the original link to gene have one</li>
    <li>The link of pathway ID. Same requirement as link to gene ID</li>
    <li>Example of gene IDs. One ID for each line. These gene IDs should be exist in the pathway annotation file.</li>
	</ul>
<li>It's recommended that you also install the <a href="http://drupal.org/project/tripal_daemon">Tripal Job Daemon</a> to manage jobs and ensure they are run soon after being submitted by the user. Without this additional module, administrators will have to execute the tripal jobs either manually or through use of cron jobs.</li>

<li>view <a href="<?php print url('pwyenrich');?>" >pathway enrichment</a> to submit a job.</li>

</ol>

