<?php
/**
 * Display the results of a BLAST job execution
 */
?>

<?php 
  // output pathway enrichment table
  if ($result_table) {

    $temp = array_keys($result_table[0]);
    array_shift($temp);
    $header_data = $temp;

    $rows_data = array(); 
    foreach ($result_table as $line) 
    {
      $line_html = array();
      $pwy_id = $line['Pathway ID'];
      $pwy_name = $line['Pathway Name'];
      

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
