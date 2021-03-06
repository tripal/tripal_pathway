<?php
/**
 * Display the results of pathway enrichment
 */

// dispaly breadcrumb
$nav = array();
$nav[] = l('Home', '/');
$nav[] = l('Pathway enrichment', 'pwyenrich');
$nav[] = t('Pathway enrichment result');
$breadcrumb = '<nav class="nav">' . implode(" > ", $nav) . '</nav><br>';
print $breadcrumb;


  // output pathway enrichment table
  if ($result_table) {

    // download link
    // dpm($tripal_args);
    $pwy_file = '../../' . $tripal_args['output_file'];
    $download_link = l('Tab-delimit format', $pwy_file );
    print '<p>Download pathway enrichment result: ' . $download_link . '</p>';

    // result table
    $temp = array_keys($result_table[0]);
    array_shift($temp);
    $header_data = $temp;

    $rows_data = array(); 
    foreach ($result_table as $line) 
    {
      $line_html = array();
      $pwy_id = $line['Pathway ID'];
      $pwy_name = $line['Pathway name'];

      $line_html[0] = "<a href=http://$dblink" . "$pwy_id target=_blank>$pwy_id<br>$pwy_name</a>"; 
      $line_html[1] = $line[$header_data[1]];
      $line_html[2] = $line[$header_data[2]];
      $line_html[3] = $line[$header_data[3]];
      
      // result with q value
      $genes = '';
      if (!empty($header_data[5])) {
         $line_html[4] = $line[$header_data[4]];
         $genes = explode(', ', $line[$header_data[5]]);
         $genes_html = '';
         foreach ($genes as $gid) {
           $genes_html.= "<a href=http://$idlink" . "$gid target=_blank>$gid</a>, ";
         }
         $line_html[5] = $genes_html;
      } 
      // result without q value
      else {
         $genes = explode(', ', $line[$header_data[4]]);
         $genes_html = '';
         foreach ($genes as $gid) {
           $genes_html.= "<a href=http://$idlink" . "$gid target=_blank>$gid</a>, ";
         }
         $line_html[4] = $genes_html;
      }

      $rows_data[] = $line_html;
    }
 
    $header = array(
      'data' => $header_data,
    );

    $rows = array(
      'data' => $rows_data,
    );


    $variables = array(
      'header' => $header_data,
      'rows' => $rows_data,
      'attributes' => array('class' => 'table'),
    );

    $table_html = theme('table', $variables);
    print $table_html;
  }
